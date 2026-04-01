local args = {}
args.size, args.constraints, args.factories = string.match(..., "^pk(%d+)(%w%w)f(%d+)$")
assert(args.size, "usage: pk{size}{sy|rc}f{fac}")

local CDA_FILE = string.format("pk%s.cda", args.size)
local XDA_FILE = string.format("pk%s.xda", args.size)
local DISTANCE_FILE = string.format("pk%s_tiem_ka_kp.txt", args.size)
local XDA_COLS_FILE = "variables.txt"
local FACTORIES_FILE = "factories_prices.txt"
local FACTORY_DISTANCE_FILE = "tiem_kp_kp.txt"
local CONSTRAINTS_FILE = string.format("con-%s.txt", args.constraints)
local MPS_FILE = (...) .. ".mps"
local OBJECTIVE_FILE = nil
local DISCOUNT = 1.04
local PERIODS = 5
local INCLUDE_CHIPS = args.factories == "72"

local df = {}
for p=1, PERIODS+1 do
	df[p] = DISCOUNT^(10*(1/2-p))
end
local dff = DISCOUNT^(-PERIODS*10)

local FIXED_OBJECTIVE = {
	NUt4 = dff,
	["WC#1"]=-df[1], ["WC#2"]=-df[2], ["WC#3"]=-df[3], ["WC#4"]=-df[4], ["WC#5"]=-df[5],
	["WEc#1"]=-df[1], ["WEc#2"]=-df[2], ["WEc#3"]=-df[3], ["WEc#4"]=-df[4], ["WEc#5"]=-df[5],
}

local CDA_COLS = {
	"ns", "KUVIONUMERO"
}
local XDA_COLS = {}
for v in io.open(XDA_COLS_FILE, "r"):lines() do
	table.insert(XDA_COLS, v)
end

local FACTORY_NODES = {
	saw1 = { type = "sawmill", capacity = "scap1" },
	saw2 = { type = "sawmill", capacity = "scap2" },
	saw3 = { type = "sawmill", capacity = "scap3" },
	pulp1 = { type = "pulpmill", capacity = "pcap1" },
	pulp2 = { type = "pulpmill", capacity = "pcap2" },
	pulp3 = { type = "pulpmill", capacity = "pcap3" },
}
if INCLUDE_CHIPS then
	FACTORY_NODES.chips_meha = { type = "chips", capacity = "meha_kim3" }
	FACTORY_NODES.chips_teolha = { type = "chips", capacity = "teolha_kim3" }
	FACTORY_NODES.homeuse = { type = "homeuse", capacity = "captot" }
end

local UNIT_FACTORY_EDGES = {
	{ to = "saw1", item = "RL1", price = "sprice1" },
	{ to = "saw2", item = "RL2", price = "sprice2" },
	{ to = "saw3", item = "RL3", price = "sprice3" },
	{ to = "pulp1", item = "RP1", price = "pprice1" },
	{ to = "pulp2", item = "RP2", price = "pprice2" },
	{ to = "pulp3", item = "RP3", price = "pprice3" },
}
if INCLUDE_CHIPS then
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "RP1", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "RP2", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "RP3", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "REst", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "REcr", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "REsr", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "chips_meha", item = "RCLP", price = "cprice2" })
	table.insert(UNIT_FACTORY_EDGES, { to = "homeuse", item = "RP3", price = "price", weight = 1/1.25, min_weight = 1/0.75 })
	table.insert(UNIT_FACTORY_EDGES, { to = "homeuse", item = "RCLP", price = "price", weight = 1/1.25, min_weight = 1/0.75 })
end

local SBPP = 0.3
local HBPP = 0.315

local FACTORY_FACTORY_EDGES = {
	{ from = "saw1", to = "pulp1", price = "pprice1" },
	{ from = "saw2", to = "pulp2", price = "pprice2" },
	{ from = "saw3", to = "pulp3", price = "pprice3" },
}
if INCLUDE_CHIPS then
	table.insert(FACTORY_FACTORY_EDGES, { from = "saw1", to = "chips_teolha", price = "cprice2" })
	table.insert(FACTORY_FACTORY_EDGES, { from = "saw2", to = "chips_teolha", price = "cprice2" })
	table.insert(FACTORY_FACTORY_EDGES, { from = "saw3", to = "chips_teolha", price = "cprice2" })
end

local FACTORY_EDGE_CONSTRAINTS = {
	{
		from = { saw1=1, saw2=1 },
		items = { RP1=1, RP2=1 },
		to = "pulpmill",
		capacity = "pcapsoftw"
	}
}

local COSTKM = {
	sawmill = 0.07,
	pulpmill = 0.08,
	chips = 0.10,
	homeuse = 0
}

--------------------------------------------------------------------------------

local ffi = require "ffi"
require "table.clear"
require "table.new"

ffi.cdef [[
	void *malloc(size_t);
	void *realloc(void *, size_t);
	void free(void *);
]]

local sizeof_u32 = 4
local sizeof_u8 = 1
local sizeof_double = 8

---- Problem definition --------------------------------------------------------

local function problem_info(_, fmt, ...)
	print(string.format("[*] "..fmt, ...))
end

local zvar_mt = { what = "zvar" }
zvar_mt.__index = zvar_mt

local function problem_add_var(problem)
	local zvar = setmetatable({z=problem.Z}, zvar_mt)
	problem.zvars[problem.Z] = zvar
	problem.Z = problem.Z + 1
	return zvar
end

local function problem_add_constraint(problem, lhs, operator, rhs)
	table.insert(problem.constraints, {lhs=lhs, operator=operator, rhs=rhs})
end

local function problem_add_unit(problem, J)
	assert(not problem.Jx, "cannot add units after adding schedules")
	problem.J[problem.I] = J
	problem.JJ[problem.I] = problem.Jtot1
	problem.Jtot1 = problem.Jtot1 + J
	problem.I = problem.I + 1
end

local function problem_cmpjbuf(problem, j1, j2)
	for k=1, problem.K do -- skip objective at 0
		local x1 = problem.jbuf[j1*(problem.K+1)+k]
		local x2 = problem.jbuf[j2*(problem.K+1)+k]
		if x1 ~= x2 then
			return x1 < x2 and "<" or ">"
		end
	end
	return "="
end

local function problem_init_schedules(problem)
	-- allocate schedule data. this overallocates because of possible duplicate schedules,
	-- but it's good enough. it's just virtual memory anyway.
	problem.Jx = ffi.gc(ffi.cast("uint32_t *", ffi.C.malloc(sizeof_u32*(problem.Jtot1+1))), ffi.C.free)
	problem.Jo = ffi.gc(ffi.cast("double *", ffi.C.malloc(sizeof_double*problem.Jtot1)), ffi.C.free)
	-- assign item numbers to items that appear in constraints
	for _,c in ipairs(problem.constraints) do
		for var in pairs(c.lhs) do
			if var.what == "item" and not var.k then
				var.k, problem.K = problem.K, problem.K+1
			end
		end
	end
	-- assign item numbers to items that appear in transports
	for t=0, problem.T-1 do
		local item = problem.transports[t].item
		if not item.k then
			item.k, problem.K = problem.K, problem.K+1
		end
	end
	problem_info(problem, "effective items: %d", problem.K)
	-- collect objective-only items
	problem.jobj = {}
	for v,c in pairs(problem.objective) do
		if v.what == "item" and not v.k then
			problem.jobj[v.name] = c
		end
	end
	-- all items that don't yet have a transport get one now
	problem.k2t = table.new(problem.K, 0)
	problem.k2name = table.new(problem.K, 0)
	for k=0, problem.K-1 do
		problem.k2t[k] = {}
	end
	for t=0, problem.T-1 do
		table.insert(problem.k2t[problem.transports[t].item.k], t)
	end
	for name,item in pairs(problem.items) do
		if item.k then
			problem.k2name[item.k] = name
			if not problem.k2t[item.k][1] then
				problem.transports[problem.T] = { item=item }
				problem.k2t[item.k][1] = problem.T
				problem.T = problem.T + 1
			end
		end
	end
	-- create utility table
	problem.Iu = ffi.gc(ffi.cast("double *", ffi.C.malloc(sizeof_double * problem.T * problem.I)),
		ffi.C.free)
	for i=0, problem.I-1 do
		for t=0, problem.T-1 do
			local util = problem.transports[t].util
			local u
			if util then
				u = util[i+1]
			else
				u = 0
			end
			problem.Iu[i*problem.T+t] = u
		end
	end
	for t=0, problem.T-1 do
		problem.transports[t].util = nil
	end
	-- set up schedule builder state
	problem.ip = 0 -- current unit
	problem.jp = 0 -- current schedule in unit
	local jjmax = 0
	for i=0, problem.I-1 do jjmax = math.max(jjmax, problem.JJ[i]) end
	problem.jbuf = ffi.gc(ffi.cast("double *", ffi.C.malloc(sizeof_double*(jjmax+1))), ffi.C.free)
	problem.jorder = {}
	-- callback for table.sort
	problem.cmpjbuf = function(j1, j2)
		local c = problem_cmpjbuf(problem, j1, j2)
		if c ~= "=" then return c == "<" end
		return problem.jbuf[j1*(problem.K+1)] > problem.jbuf[j2*(problem.K+1)]
	end
end

local function problem_grownz(problem)
	ffi.gc(problem.Xk, nil)
	ffi.gc(problem.Xx, nil)
	problem.sizenz = math.floor(1.5 * math.max(problem.sizenz, problem.Jtot1))
	problem_info(problem, "grow nonzero array (%d million entries = %dMB), read %d units",
		problem.sizenz/10^6, problem.sizenz*(sizeof_u8+sizeof_double)/2^20, problem.ip)
	problem.Xk = ffi.gc(ffi.cast("uint8_t *", ffi.C.realloc(problem.Xk, problem.sizenz*sizeof_u8)),
		ffi.C.free)
	problem.Xx = ffi.gc(ffi.cast("double *", ffi.C.realloc(problem.Xx, problem.sizenz*sizeof_double)),
		ffi.C.free)
end

local function problem_commit_schedules(problem)
	table.clear(problem.jorder)
	for i=1, problem.J[problem.ip] do
		problem.jorder[i] = i-1
	end
	table.sort(problem.jorder, problem.cmpjbuf)
	problem.JJ[problem.ip] = problem.Jtot
	for jj=1, problem.J[problem.ip] do
		local j = problem.jorder[jj]
		if jj == 1 or problem_cmpjbuf(problem, j, problem.jorder[jj-1]) ~= "=" then
			problem.Jo[problem.Jtot] = problem.jbuf[j*(problem.K+1)] -- objective
			problem.Jx[problem.Jtot] = problem.nnz -- nonzeros start
			problem.Jtot = problem.Jtot+1
			for k=0, problem.K-1 do
				local x = problem.jbuf[j*(problem.K+1)+k+1]
				if x ~= 0 then
					if problem.nnz >= problem.sizenz then
						problem_grownz(problem)
					end
					problem.Xk[problem.nnz] = k
					problem.Xx[problem.nnz] = x
					problem.nnz = problem.nnz + 1
				end
			end
		end
	end
	problem.J[problem.ip] = problem.Jtot - problem.JJ[problem.ip]
	problem.ip = problem.ip+1
	problem.jp = 0
end

local function problem_add_schedule(problem, row)
	if not problem.Jx then
		problem_init_schedules(problem)
	end
	local obj = 0
	for v,c in pairs(problem.jobj) do
		obj = obj + c*row[v]
	end
	problem.jbuf[problem.jp*(problem.K+1)] = obj
	for k=0, problem.K-1 do
		problem.jbuf[problem.jp*(problem.K+1)+k+1] = row[problem.k2name[k]]
	end
	problem.jp = problem.jp+1
	if problem.jp == problem.J[problem.ip] then
		problem_commit_schedules(problem)
	end
end

local function problem_finishschedules(problem)
	if problem.ip then
		if problem.ip ~= problem.I then
			error(string.format("expected %d units, got %d", problem.I, problem.ip))
		end
		problem.Jx[problem.Jtot] = problem.nnz
		problem_info(problem, "%d nonzeros (%d allocated) = %dMB (%dMB)", problem.nnz,
			problem.sizenz, problem.nnz*(sizeof_u8+sizeof_double)/2^20,
			problem.sizenz*(sizeof_u8+sizeof_double)/2^20)
		problem_info(problem, "sparsity: %g", 1-problem.nnz/(problem.Jtot*problem.K))
		problem_info(problem, "effective schedules: %d (removed %d)", problem.Jtot,
			problem.Jtot1-problem.Jtot)
		problem.jbuf = nil
		problem.jorder = nil
		problem.ip = nil
		problem.jp = nil
	end
end

local item_mt = { what = "item" }
item_mt.__index = item_mt

local function problem_item(problem, item)
	local k = problem.items[item]
	if not k then
		assert(not problem.Jx, "cannot add items after adding schedules")
		k = setmetatable({name=item}, item_mt)
		problem.items[item] = k
	end
	return k
end

local transport_mt = { what = "transport" }
transport_mt.__index = transport_mt

local function problem_add_transport(problem, item, util)
	assert(not problem.Jx, "cannot add transports after adding schedules")
	if #util ~= problem.I then
		error(string.format("expected %d utilities, got %d", problem.I, #util))
	end
	local transport = {
		item = problem_item(problem, item),
		util = util,
		what = "transport",
		t = problem.T
	}
	problem.transports[problem.T] = transport
	problem.T = problem.T + 1
	return transport
end

---- MPS Writer ----------------------------------------------------------------

-- rows:
--   W{i}          → unit i area constraint                  ∑_j w_ij = 1
--   X{i}.{k}      → unit i, item k transport constraint     ∑_j r_ijk*w_ij - ∑_t x_it = 0
--   T{t}          → transport t definition (implicit k)     t_t - ∑_i x_it = 0
--   K{k}          → item k definition                       k_k - ∑_t t_kt = 0
--   R{i}          → auxiliary row i
-- columns:
--   w{i}.{j}      → unit i schedule j
--   x{i}.{t}      → unit i, transport t (item k implicit)
--   {item}%{t}    → item `item`, transport t (total)
--   {item}        → item `item` (total)
--   z{i}          → zvar i
-- bounds:
--   w{i}.{j}      ≥ 0
--   x{i}.{t}      ≥ 0   if item k has at least 2 transports
--   x{i}.{t}      free  if item k has exactly one transport
--   {item}%{t}    free
--   {item}        free
--   z{i}          ≥ 0

local mps_rowtype = { ["="] = "E", ["<"] = "L", [">"] = "G" }

local function mps_rows(problem, fp)
	-- area constraints Wi:
	for i=0, problem.I-1 do
		fp:write(" E W", i, "\n")
	end
	-- transport constraints Xi.k:
	for k=0, problem.K-1 do
		for i=0, problem.I-1 do
			fp:write(" E X", i, ".", k, "\n")
		end
	end
	-- transport totals Tt:
	for k=0, problem.K-1 do
		for _,t in ipairs(problem.k2t[k]) do
			fp:write(" E T", t, "\n")
		end
	end
	-- item totals Kk:
	for k=0, problem.K-1 do
		fp:write(" E K", k, "\n")
	end
	-- auxiliary rows Ri:
	for name,row in pairs(problem.mps.rows) do
		fp:write(" ", row.op, " ", name, "\n")
	end
end

local function mps_put(mps, fp, value, name1, ...)
	if value == 0 then
		return
	end
	if mps.cur == name1 then
		mps.cur = nil
	else
		mps.cur = name1
		fp:write("\n ", name1)
	end
	fp:write(" ")
	fp:write(...)
	fp:write(" ", value)
end

local function mps_columns(problem, fp)
	-- weights wi.j
	for i=0, problem.I-1 do
		for j=0, problem.J[i]-1 do
			local jj = problem.JJ[i]+j
			local wij = string.format("w%d.%d", i, j)
			mps_put(problem.mps, fp, 1, wij, "W", i)
			mps_put(problem.mps, fp, problem.Jo[jj]*(1/problem.scale), wij, "obj")
			for v=problem.Jx[jj], problem.Jx[jj+1]-1 do
				mps_put(problem.mps, fp, problem.Xx[v]*(1/problem.scale), wij, "X", i, ".", problem.Xk[v])
			end
		end
	end
	-- transports xi.t
	for k=0, problem.K-1 do
		for _,t in ipairs(problem.k2t[k]) do
			for i=0, problem.I-1 do
				local xit = string.format("x%d.%d", i, t)
				-- local tx = problem.transports[t]
				-- if tx.util then
				-- 	mps_put(problem.mps, fp, tx.util[i+1], xit, "obj")
				-- end
				mps_put(problem.mps, fp, problem.Iu[i*problem.T+t], xit, "obj")
				mps_put(problem.mps, fp, -1, xit, "T", t)
				mps_put(problem.mps, fp, -1, xit, "X", i, ".", k)
			end
		end
	end
	-- item transports {item}%{t}, total transports {item}, and zvars zi
	for colname,col in pairs(problem.mps.cols) do
		for rowname,value in pairs(col) do
			mps_put(problem.mps, fp, value, colname, rowname)
		end
	end
	fp:write("\n")
end

local function mps_rhs(problem, fp)
	-- area constraints Wi
	for i=0, problem.I-1 do
		mps_put(problem.mps, fp, problem.scale, "RHS", "W", i)
	end
	-- Xik rhs=0
	-- Tt rhs=0
	-- Kk rhs=0
	-- auxiliary rows Ri
	for name,row in pairs(problem.mps.rows) do
		mps_put(problem.mps, fp, row.rhs, "RHS", name)
	end
	fp:write("\n")
end

local function mps_bounds(problem, fp)
	-- wi.j ≥ 0
	-- xi.t free if #k2t[k] = 1
	for k=0, problem.K-1 do
		if problem.k2t[k][2] == nil then
			local t = problem.k2t[k][1]
			for i=0, problem.I-1 do
				fp:write(" FR BND x", i, ".", t, "\n")
			end
		end
	end
	-- {item}%{t}, {item} free, zi ≥ 0
	for name in pairs(problem.mps.cols) do
		if not name:find("^z%d+$") then
			fp:write(" FR BND ", name, "\n")
		end
	end
end

local function mps_col(mps, var)
	if type(var) == "string" then
		return var
	elseif var.what == "item" then
		return var.name
	elseif var.what == "transport" then
		return string.format("%s%%%d", var.item.name, var.t)
	elseif var.what == "zvar" then
		local zv = string.format("z%d", var.z)
		if not mps.cols[zv] then
			mps.cols[zv] = {}
		end
		return zv
	else
		error("TODO")
	end
end

local function mps_init(problem)
	problem.mps = {
		rows = {}, -- idx → { op="E"|"L"|"G", rhs=number }
		cols = {}, -- name → { [row]=value }
		items = {}
	}
	for k=0, problem.K-1 do
		local name = problem.k2name[k]
		problem.mps.cols[name] = { [string.format("K%d", k)]=1 }
		for _,t in ipairs(problem.k2t[k]) do
			problem.mps.cols[string.format("%s%%%d", name, t)] = {
				[string.format("K%d", k)] = -1,
				[string.format("T%d", t)] = 1
			}
		end
	end
	for i,c in ipairs(problem.constraints) do
		local row = string.format("R%d", i)
		problem.mps.rows[row] = { op=mps_rowtype[c.operator], rhs=c.rhs }
		for var,coef in pairs(c.lhs) do
			problem.mps.cols[mps_col(problem.mps, var)][row] = coef
		end
	end
	for var,coef in pairs(problem.objective) do
		if var.what ~= "item" or var.k then
			problem.mps.cols[mps_col(problem.mps, var)].obj = coef
		end
	end
end

local function problem_write_mps(problem, fp, name)
	problem.scale = problem.scale or 1
	problem_finishschedules(problem)
	mps_init(problem)
	fp:write("NAME ", name or "fprob", "\n")
	fp:write("OBJSENSE MAX\n")
	fp:write("ROWS\n")
	fp:write(" N obj\n")
	mps_rows(problem, fp)
	fp:write("COLUMNS") -- no \n
	mps_columns(problem, fp)
	fp:write("RHS") -- no \n
	mps_rhs(problem, fp)
	fp:write("BOUNDS\n")
	mps_bounds(problem, fp)
	fp:write("ENDATA\n")
	problem.mps = nil
end

--------------------------------------------------------------------------------

local problem_mt = {
	info            = problem_info,
	add_var         = problem_add_var,
	add_constraint  = problem_add_constraint,
	add_unit        = problem_add_unit,
	add_schedule    = problem_add_schedule,
	add_transport   = problem_add_transport,
	item            = problem_item,
	write_mps       = problem_write_mps,
}
problem_mt.__index = problem_mt

local function newproblem()
	return setmetatable({
		I             = 0,     -- number of units
		Iu            = nil,   -- unit*T + transport -> util
		J             = {},    -- unit -> number of schedules
		JJ            = {},    -- unit -> index of first schedule (0-based)
		Jtot1         = 0,     -- total number of schedules
		Jtot          = 0,     -- number of effective schedules
		K             = 0,     -- number of items
		T             = 0,     -- number of transports
		Z             = 0,     -- number of zvars
		items         = {},    -- item name -> {k, name}
		transports    = {},    -- transport number -> {item, util, t}
		objective     = {},    -- objective row
		constraints   = {},    -- auxiliary constraints
		Jx            = nil,   -- schedule number -> start of nonzeros   (uint32_t)
		Jo            = nil, -- schedule number -> objective  (double)
		Xk            = ffi.cast("uint8_t *", 0),   -- nonzero number -> item number
		Xx            = ffi.cast("double *", 0),   -- nonzero number -> value
		nnz           = 0,     -- number of nonzeros
		sizenz        = 0,     -- size of nonzero array
		zvars         = {},    -- zvar number -> {z}
	}, problem_mt)
end

--------------------------------------------------------------------------------

ffi.cdef [[
	double strtod(const char *, char **);
]]

local function openread(fp)
	if type(fp) == "string" then
		return assert(io.open(fp, "r"))
	else
		return fp
	end
end

local strep = ffi.new("char *[1]")
local function read_matrix_row(ctx)
	local line = ctx.fp:read("*l")
	if not line then
		ctx.fp:close()
		return
	end
	local len = #line
	local ptr = ffi.cast("const char *", line)
	local ofs = line:find("[^%s]")-1
	local row = ctx.row
	if row then
		table.clear(row)
	else
		row = {}
	end
	for i=1, ctx.ncol do
		row[ctx.cols[i]] = ffi.C.strtod(ptr+ofs, strep)
		local e = strep[0] - ptr + 1
		if e == len then
			assert(i == ctx.ncol, "too few columns")
			break
		end
		ofs = e
		-- local s,e = line:find("%s+", ofs)
		-- -- row[ctx.cols[i]] = tonumber(line:sub(ofs,s))
		-- if not e then
		-- 	assert(i == ctx.ncol, "too few columns")
		-- 	break
		-- end
		-- ofs = e+1
	end
	return row
end

local function read_matrix(fp, cols, rowbuf)
	return read_matrix_row, {fp=openread(fp), cols=cols, ncol=#cols, spaces={}, row=rowbuf}
end

--------------------------------------------------------------------------------

local jfiles = require "jfiles"
require "table.new"

local var2name = {}
local problem = newproblem()
for item,c in pairs(FIXED_OBJECTIVE) do
	local var = problem:item(item)
	var2name[var] = item
	problem.objective[var] = c
end

print(string.format("[*] reading units [%s]", CDA_FILE))
local units = {}
do
	local ns = 0
	for row in read_matrix(CDA_FILE, CDA_COLS) do
		problem:add_unit(row.ns)
		ns = ns + row.ns
		table.insert(units, row)
	end
	print(string.format("--> %d units, %d schedules", #units, ns))
end

if FACTORIES_FILE then
	print(string.format("[*] reading factories [%s]", FACTORIES_FILE))
	local factories = jfiles.read_factories(FACTORIES_FILE)
	factories.homeuse[1].price = 15
	print(string.format("[*] reading distance matrix [%s]", DISTANCE_FILE))
	local distance = jfiles.read_distance(DISTANCE_FILE)
	for _,u in ipairs(units) do distance[u.KUVIONUMERO].ppESKT = 0 end
	local factorydistance
	if FACTORY_DISTANCE_FILE then
		print(string.format("[*] reading factory-factory distance matrix [%s]", FACTORY_DISTANCE_FILE))
		factorydistance = jfiles.read_distance(FACTORY_DISTANCE_FILE, true)
	end
	print("[*] computing transports")
	local nnode, nufe, nffe, nfc, nxc = 0, 0, 0, 0, 0
	local nodes = {}
	for _,fs in pairs(factories) do
		for _,f in ipairs(fs) do
			f.nodes = {}
			for p=1, PERIODS do
				f.nodes[p] = {}
			end
		end
	end
	for ntype, info in pairs(FACTORY_NODES) do
		nodes[ntype] = {}
		for p=1, PERIODS do
			nodes[ntype][p] = {}
		end
		for _,factory in ipairs(factories[info.type]) do
			local cap = type(info.capacity) == "number" and info.capacity or factory[info.capacity]
			if cap > 0 then
				for p=1, PERIODS do
					local node = {factory=factory, inputs={}, outputs={}}
					table.insert(nodes[ntype][p], node)
					table.insert(factory.nodes[p], node)
					nnode = nnode+1
				end
			end
		end
	end
	print(string.format("--> %d factory nodes", nnode))
	for _,edge in ipairs(UNIT_FACTORY_EDGES) do
		for p=1, PERIODS do
			for _,node in ipairs(nodes[edge.to][p]) do
				local price = node.factory[edge.price]
				local costkm = COSTKM[FACTORY_NODES[edge.to].type]
				local util = table.new(#units, 0)
				for i=1, #units do
					local dist = distance[units[i].KUVIONUMERO][node.factory.name]
					util[i] = (price - costkm*dist) * df[p]
				end
				local var = problem:add_transport(string.format("%s#%d", edge.item, p), util)
				var2name[var] = string.format("%s_%d%%%%%s", edge.item, p, node.factory.name)
				table.insert(node.inputs, {edge=edge, var=var})
				nufe = nufe+1
			end
		end
	end
	print(string.format("--> %d unit-factory edges", nufe))
	if factorydistance then
		for _,edge in ipairs(FACTORY_FACTORY_EDGES) do
			for p=1, PERIODS do
				for _,sink in ipairs(nodes[edge.to][p]) do
					for _,source in ipairs(nodes[edge.from][p]) do
						local price = sink.factory[edge.price]
						local costkm = COSTKM[FACTORY_NODES[edge.to].type]
						local dist = factorydistance[source.factory.name][sink.factory.name]
						local var = problem:add_var()
						var2name[var] = string.format("RL%d_%d%%%s__%s", edge.from:sub(4,4), p,
						source.factory.name, sink.factory.name)
						problem.objective[var] = (price - costkm*dist) * df[p]
						local e = {edge=edge, var=var}
						table.insert(source.outputs, e)
						table.insert(sink.inputs, e)
						nffe = nffe+1
						--source.bp_source = true
						--sink.bp_sink = true
					end
				end
			end
		end
		print(string.format("--> %d factory-factory edges", nffe))
	end
	for ntype,nx in pairs(nodes) do
		for p=1, PERIODS do
			for _,node in ipairs(nx[p]) do
				if node.inputs[1] then
					local row = {}
					local lowbound
					for _,e in ipairs(node.inputs) do
						row[e.var] = e.edge.weight or 1
						if e.edge.min_weight then
							if not lowbound then lowbound = {} end
							lowbound[e.var] = e.edge.min_weight
						end
					end
					local cap = FACTORY_NODES[ntype].capacity
					if type(cap) == "string" then
						cap = node.factory[cap]
					end
					if cap < math.huge and not node.bp_sink then
						problem:add_constraint(row, "<", cap)
						--for v,c in pairs(row) do
						--	io.stdout:write(c, "*", var2name[v], "\n")
						--end
						--io.stdout:write("< ", node.factory.name, "%", FACTORY_NODES[ntype].capacity, "\n\n")
						nfc = nfc+1
					end
					if lowbound and cap < math.huge then
						problem:add_constraint(lowbound, ">", cap)
						nfc = nfc+1
					end
					if node.outputs[1] and not node.bp_source then
						local row = {}
						for _,e in ipairs(node.inputs) do
							-- row[e.var] = e.edge.weight or 1
							row[e.var] = var2name[e.var]:sub(1,3) == "RL3" and HBPP or SBPP
						end
						for _,e in ipairs(node.outputs) do
							row[e.var] = -(e.edge.weight or 1)
						end
						problem:add_constraint(row, ">", 0)
						--for v,c in pairs(row) do
						--	io.stdout:write(c, "*", var2name[v], "\n")
						--end
						--io.stdout:write("> 0\n\n")
						nfc = nfc+1
					end
				end
			end
		end
	end
	print(string.format("--> %d flow constraints", nfc))
	for _,con in ipairs(FACTORY_EDGE_CONSTRAINTS) do
		for _,f in ipairs(factories[con.to]) do
			for p=1, PERIODS do
				local lhs = {}
				for _,node in ipairs(f.nodes[p]) do
					for _,input in ipairs(node.inputs) do
						if input.edge.item and con.items[input.edge.item] then
							lhs[input.var] = con.items[input.edge.item]
						elseif input.edge.from and con.from[input.edge.from] then
							lhs[input.var] = con.from[input.edge.from]
						end
					end
				end
				if next(lhs) then
					--for v,c in pairs(lhs) do
					--	io.stdout:write(c, "*", var2name[v], "\n")
					--end
					--io.stdout:write("< ", f.name, "%", con.capacity,"\n\n")
					problem:add_constraint(lhs, "<", f[con.capacity])
					nxc = nxc+1
				end
			end
		end
	end
	print(string.format("--> %d factory-edge constraints", nxc))
end

if OBJECTIVE_FILE then
	print(string.format("[*] reading objective [%s]", OBJECTIVE_FILE))
	for _,term in ipairs(jfiles.read_expr(OBJECTIVE_FILE)) do
		local var = problem:item(term.var)
		var2name[var] = term.var
		problem.objective[var] = term.coef or 1
	end
end

if CONSTRAINTS_FILE then
	print(string.format("[*] reading constraints [%s]", CONSTRAINTS_FILE))
	for _,row in ipairs(jfiles.read_constraints(CONSTRAINTS_FILE)) do
		local lhs = {}
		for _,term in ipairs(row.lhs) do
			local var = problem:item(term.var)
			var2name[var] = term.var
			lhs[var] = term.coef or 1
		end
		assert(not row.rhs.var, "TODO?")
		problem:add_constraint(lhs, row.operator, row.rhs.coef)
	end
end

print(string.format("[*] reading schedules [%s]", XDA_FILE))
for row in read_matrix(XDA_FILE, XDA_COLS, table.new(0, #XDA_COLS)) do
	problem:add_schedule(row)
end

--------------------------------------------------------------------------------

problem.scale = tonumber(args.size)
print(string.format("[*] writing mps [%s]", MPS_FILE))
local fp = io.open(MPS_FILE, "w")
problem:write_mps(fp)
fp:close()
