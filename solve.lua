local args = {}
args.size, args.constraints, args.factories, args.mode = string.match(..., "^pk(%d+)(%w%w)f(%d+)(%a)$")
assert(args.size, "usage: pk{size}{sy|rc}f{fac}{mode}")

local CDA_FILE = string.format("pk%s.cda", args.size)
local XDA_FILE = string.format("pk%s.xda", args.size)
local DISTANCE_FILE = string.format("pk%s_tiem_ka_kp.txt", args.size)
local XDA_COLS_FILE = "variables.txt"
local FACTORIES_FILE = "factories_prices.txt"
local FACTORY_DISTANCE_FILE = "tiem_kp_kp.txt"
local CONSTRAINTS_FILE = string.format("con-%s.txt", args.constraints)
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

math.randomseed(123)

--------------------------------------------------------------------------------

local jfiles = require "jfiles"
require "table.new"

local distance
local function util(unit, transport)
	if not distance then return 0 end
	return df[transport.period]
		* (transport.price - transport.costkm*distance[unit.KUVIONUMERO][transport.factory])
end

local problem = {
	threads     = 12,
	mode        = args.mode,
	cda         = CDA_FILE,
	cda_cols    = CDA_COLS,
	items       = {},
	xda         = XDA_FILE,
	transports  = {},
	constraints = {},
	objective   = FIXED_OBJECTIVE,
	utility     = util,
	info        = function(x) io.write("[*] ", x, "\n") end,
}
for v in io.open(XDA_COLS_FILE, "r"):lines() do
	table.insert(problem.items, v)
end

if FACTORIES_FILE then
	print(string.format("[*] reading factories [%s]", FACTORIES_FILE))
	local factories = jfiles.read_factories(FACTORIES_FILE)
	factories.homeuse[1].price = 15
	print(string.format("[*] reading distance matrix [%s]", DISTANCE_FILE))
	distance = jfiles.read_distance(DISTANCE_FILE)
	for _,r in pairs(distance) do r.ppESKT = 0 end
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
				local item = string.format("%s#%d", edge.item, p)
				local factory = node.factory.name
				local tx = {
					item = item,
					factory = factory,
					costkm = COSTKM[FACTORY_NODES[edge.to].type],
					name = string.format("%s%%%%%s", item, factory),
					price = node.factory[edge.price],
					period = p
				}
				nufe = nufe+1
				table.insert(problem.transports, tx)
				table.insert(node.inputs, {edge=edge, var=tx})
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
						local var = string.format("RL%d_%d%%%s__%s", edge.from:sub(4,4), p,
							source.factory.name, sink.factory.name)
						problem.objective[var] = (price - costkm*dist) * df[p]
						local e = {edge=edge, var=var}
						table.insert(source.outputs, e)
						table.insert(sink.inputs, e)
						nffe = nffe+1
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
						table.insert(problem.constraints, {lhs=row, sense="<", rhs=cap})
						-- for v,c in pairs(row) do
						-- 	io.stdout:write(c, "*", type(v) == "string" and v or v.name, "\n")
						-- end
						-- io.stdout:write("< ", node.factory.name, "%", FACTORY_NODES[ntype].capacity, "\n\n")
						nfc = nfc+1
					end
					if lowbound and cap < math.huge then
						table.insert(problem.constraints, {lhs=lowbound, sense=">", rhs=cap})
						nfc = nfc+1
					end
					if node.outputs[1] and not node.bp_source then
						local row = {}
						for _,e in ipairs(node.inputs) do
							-- row[e.var] = e.edge.weight or 1
							row[e.var] = e.var.name:sub(1,3) == "RL3" and HBPP or SBPP
						end
						for _,e in ipairs(node.outputs) do
							row[e.var] = -(e.edge.weight or 1)
						end
						table.insert(problem.constraints, {lhs=row, sense=">", rhs=0})
						-- for v,c in pairs(row) do
						-- 	io.stdout:write(c, "*", type(v) == "string" and v or v.name, "\n")
						-- end
						-- io.stdout:write("> 0\n\n")
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
					-- for v,c in pairs(lhs) do
					-- 	io.stdout:write(c, "*", type(v) == "string" and v or v.name, "\n")
					-- end
					-- io.stdout:write("< ", f.name, "%", con.capacity,"\n\n")
					table.insert(problem.constraints, {lhs=lhs, sense="<", rhs=f[con.capacity]})
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
		problem.objective[term.var] = term.coef or 1
	end
end

if CONSTRAINTS_FILE then
	print(string.format("[*] reading constraints [%s]", CONSTRAINTS_FILE))
	for _,row in ipairs(jfiles.read_constraints(CONSTRAINTS_FILE)) do
		local lhs = {}
		for _,term in ipairs(row.lhs) do
			lhs[term.var] = term.coef or 1
		end
		assert(not row.rhs.var, "TODO?")
		table.insert(problem.constraints, {lhs=lhs, sense=row.operator, rhs=row.rhs.coef})
	end
end

--------------------------------------------------------------------------------

require("fscg").solve(problem)
