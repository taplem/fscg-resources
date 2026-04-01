local ffi = require "ffi"
local buffer = require "string.buffer"
local CFTAB, CDEF, SOLVER = ...
ffi.cdef(CDEF)
require "table.new"

-- order must match fscg.c!
local C = ffi.cast([[
	const struct {
		void *(*malloc)(size_t);
		FSCGFile *(*fscg_file_open)(const char *);
		const char *(*fscg_file_line)(FSCGFile *, size_t *);
		int (*fscg_file_read)(FSCGFile *, int, double *);
		void (*fscg_file_close)(FSCGFile *);
		bool (*fscg_model_read_xda)(FSCGModel *, const char *, int, int *, double *);
		bool (*fscg_model_solve)(FSCGModel *);
		void (*fscg_model_destroy)(FSCGModel *);
		int (*fscg_cpus)(void);
		double (*fscg_time)(void);
	} *
]], CFTAB)

local function nextrow(ctx)
	local n = C.fscg_file_read(ctx.fp, ctx.ncol, ctx.buf)
	if n == 0 then
		return
	end
	if n < ctx.ncol then
		error(string.format("not enough columns (%d < %d)", n, ctx.ncol))
	end
	local row = table.new(0, ctx.ncol)
	for i=1, ctx.ncol do
		row[ctx.cols[i]] = ctx.buf[i-1]
	end
	return row
end

local function rows(fname, cols)
	local fp = C.fscg_file_open(fname)
	if fp == nil then
		error(string.format("failed to read `%s'", fname))
	end
	fp = ffi.gc(fp, C.fscg_file_close)
	if not cols then
		local len = ffi.new("size_t[1]")
		local p = C.fscg_file_line(fp, len)
		local header = ffi.string(p, len[0])
		cols = {}
		for name in header:gmatch("([%w%.%-_]+)") do
			table.insert(cols, name)
		end
	end
	local ctx = {
		fp   = fp,
		ncol = #cols,
		cols = cols,
		buf  = ffi.new("double[?]", #cols)
	}
	return nextrow, ctx
end

local function read(fname, cols)
	local dat = {}
	local n = 0
	for row in rows(fname, cols) do
		n = n+1
		dat[n] = row
	end
	return dat
end

local function check(m, r)
	if not r then
		error(ffi.string(m.err))
	end
end

local function noinfo() end

local function definemodel(model)
	model.info("preprocessing model")
	local m = ffi.gc(ffi.new("FSCGModel"), C.fscg_model_destroy)
	local name2col = {} -- item name -> xda column
	local name2z = {} -- z-variable name -> z
	local col2k = ffi.new("int[?]", #model.items) -- xda column -> k
	local col2o = ffi.new("double[?]", #model.items) -- xda column -> objective
	local k2item = {} -- k -> item
	-- info callback
	local model_info = model.info
	if model_info then
		model.info = function(...) return model_info(string.format(...)) end
	else
		model.info = noinfo
	end
	-- collect column numbers
	for i,name in ipairs(model.items) do
		name2col[name] = i-1
		col2k[i-1] = -1
	end
	-- collect transport items
	if model.transports then
		for _,t in ipairs(model.transports) do
			local k = col2k[name2col[t.item]]
			if k == -1 then
				k, m.K = m.K, m.K+1
				col2k[name2col[t.item]] = k
				k2item[k] = { transports={} }
			end
			table.insert(k2item[k].transports, t)
		end
	end
	-- collect objective row
	if model.objective then
		for x,v in pairs(model.objective) do
			if name2col[x] then
				col2o[name2col[x]] = v
			elseif type(x) == "string" and not name2z[x] then
				name2z[x], m.Z = m.Z, m.Z+1
			end
		end
	end
	-- collect constraint rows
	if model.constraints then
		for _,c in ipairs(model.constraints) do
			for x in pairs(c.lhs) do
				if name2col[x] then
					local k = col2k[name2col[x]]
					if k == -1 then
						k, m.K = m.K, m.K+1
						col2k[name2col[x]] = k
						k2item[k] = {}
					end
				elseif type(x) == "string" and not name2z[x] then
					name2z[x], m.Z = m.Z, m.Z+1
				end
			end
		end
	end
	-- collect transports
	local transport2t = {} -- transport object -> t
	m.Kt = C.malloc(m.K*ffi.sizeof("int"))
	m.Za = m.Z
	for k=0, m.K-1 do
		local item = k2item[k]
		item.t = m.T
		local nt
		if item.transports then
			if #item.transports >= 2 then
				item.z, m.Z = m.Z, m.Z+1
			end
			for t,tx in ipairs(item.transports) do
				transport2t[tx] = m.T+t-1
			end
			nt = #item.transports
		else
			nt = 1
		end
		m.Kt[k] = nt
		m.T = m.T+nt
	end
	-- translate objective
	local objective = {}
	local onnz = 0
	if model.objective then
		for x,v in pairs(model.objective) do
			if name2z[x] then
				objective[m.T+name2z[x]] = v
				onnz = onnz+1
			elseif transport2t[x] then
				-- TODO: move it to util
				error("TODO")
			end
		end
	end
	-- translate constraints
	local constraints = {}
	for i,c in ipairs(model.constraints) do
		local lhs = {}
		for x,v in pairs(c.lhs) do
			local k
			if name2col[x] then
				local item = k2item[col2k[name2col[x]]]
				k = item.z and (m.T+item.z) or item.t
			elseif name2z[x] then
				k = m.T+name2z[x]
			else
				k = transport2t[x]
			end
			lhs[k] = v
		end
		constraints[i] = {lhs=lhs, sense=c.sense, rhs=c.rhs}
	end
	-- define artificial z-variables
	for k=0, m.K-1 do
		local item = k2item[k]
		if item.z then
			local lhs = { [m.T+item.z]=-1 }
			for t=0, m.Kt[k]-1 do
				lhs[item.t+t] = 1
			end
			table.insert(constraints, {lhs=lhs, sense="=", rhs=0})
		end
	end
	-- define constraint matrix
	local cnnz = 0
	for _,c in ipairs(constraints) do
		for _ in pairs(c.lhs) do
			cnnz = cnnz+1
		end
	end
	m.C = #constraints
	m.CX = onnz+cnnz
	m.Ca = #model.constraints
	m.Cnnz = C.malloc((m.C+1)*ffi.sizeof("int"))
	m.Cs = C.malloc((m.C+1)*ffi.sizeof("char"))
	m.Cr = C.malloc((m.C+1)*ffi.sizeof("double"))
	m.CXc = C.malloc(m.CX*ffi.sizeof("int"))
	m.CXv = C.malloc(m.CX*ffi.sizeof("double"))
	m.Cnnz[0] = onnz
	m.Cs[0] = string.byte "!"
	m.Cr[0] = 0/0
	local nzi = 0
	for i,v in pairs(objective) do
		m.CXc[nzi] = i
		m.CXv[nzi] = v
		nzi = nzi+1
	end
	for i,c in ipairs(constraints) do
		local nzi1 = nzi
		for j,v in pairs(c.lhs) do
			m.CXc[nzi] = j
			m.CXv[nzi] = v
			nzi = nzi+1
		end
		m.Cnnz[i] = nzi-nzi1
		m.Cs[i] = string.byte(c.sense)
		m.Cr[i] = c.rhs
	end
	-- read unit data
	if not model.units then
		local cda = assert(model.cda, "model doesn't define `units' or `cda'")
		local cda_cols = assert(model.cda_cols, "model defines `cda' but not `cda_cols'")
		model.info("reading units [%s]", cda)
		model.units = read(cda, cda_cols)
	end
	-- compute utilities and schedule counts
	local txs = {}
	local txs_t = {}
	local txs_T = 0
	for tx,t in pairs(transport2t) do
		txs[txs_T] = tx
		txs_t[txs_T] = t
		txs_T = txs_T+1
	end
	m.I = #model.units
	m.Ij = C.malloc((m.I+1)*ffi.sizeof("uint32_t"))
	m.ITu = C.malloc(m.I*m.T*ffi.sizeof("double"))
	ffi.fill(m.ITu, m.I*m.T*ffi.sizeof("double"))
	local J = 0
	for i,u in ipairs(model.units) do
		local ns = assert(u.ns, "unit doesn't define `ns'")
		J = J+ns
		m.Ij[i-1] = ns
		for tt=0, txs_T-1 do
			m.ITu[(i-1)*m.T+txs_t[tt]] = model.utility(u,txs[tt])
		end
	end
	m.J = J
	model.J = J -- for summary
	-- read schedule data
	model.info("reading schedules [%s]", model.xda)
	local xda = assert(model.xda, "model doesn't define `xda'")
	m.scale = model.scale or 1
	local t = C.fscg_time()
	check(m, C.fscg_model_read_xda(m, xda, #model.items, col2k, col2o))
	model.t_read = C.fscg_time() - t
	model.info("read in %gs", model.t_read)
	return m
end

local function summary(model, m)
	local buf = buffer.new()
	buf:putf("units:           %d\n", m.I)
	buf:putf("schedules:       %d (out of %d)\n", m.J, model.J)
	buf:putf("items:           %d (out of %d)\n", m.K, #model.items)
	buf:putf("transports:      %d (%d artificial)\n", m.T, m.T-#model.transports)
	buf:putf("constraints:     %d (%d artificial)\n", m.C, m.C-#model.constraints)
	buf:putf("z-variables:     %d (%d artificial)\n", m.Z, m.Z-m.Za)
	buf:putf("model nonzeros:  %d\n", m.CX)
	buf:putf("data nonzeros:   %d\n", m.Jx[m.J])
	return buf:get()
end

local function configure(m, model)
	if model.threads == true then
		m.threads = C.fscg_cpus()
	elseif type(model.threads) == "number" then
		m.threads = model.threads
	elseif not model.threads then
		m.threads = 0
	else
		error(string.format("bad thread count: %s", model.threads))
	end
	if model.mode then
		m.mode = string.byte(model.mode)
	elseif m.threads > 1 then
		m.mode = string.byte "D"
	else
		m.mode = string.byte "p"
	end
	if model.epsilon then m.epsilon = model.epsilon end
	if model.gamma then m.gamma = model.gamma end
	if model.N0 then m.s_N0 = model.N0 end
	if model.N1 then m.s_N1 = model.N1 end
	if model.W then m.s_W = model.W end
end

-- model:
--   units = [unit]  OR  cda=string, cda_cols = [string]
--   xda = string
--   items = [string]
--   transports = [{item=..., ...}]
--   constraints = [{lhs=linexpr, sense='<'|'>'|'=', rhs=number}]
--   objective = linexpr
--   utility = function unit,transport -> number
local function solve(model)
	local m = definemodel(model)
	model.info("model summary:\n%s", summary(model, m))
	configure(m, model)
	local t = C.fscg_time()
	check(m, C.fscg_model_solve(m))
	model.t_solve = C.fscg_time() - t
	model.info("solved in %gs", model.t_solve)
end

return {
	rows   = rows,
	read   = read,
	solve  = solve,
	solver = SOLVER
}
