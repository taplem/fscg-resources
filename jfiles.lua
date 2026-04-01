local fscg = require "fscg"
local buffer = require "string.buffer"
require "table.clear"
local assert, tonumber = assert, tonumber

local function openread(fp)
	if type(fp) == "string" then
		return assert(io.open(fp, "r"))
	else
		return fp
	end
end

local function read_factories(fname)
	local fp = openread(fname)
	local factories = {}
	local state = "new"
	local cur, props
	for line in fp:lines() do
		if state == "new" then
			cur = line:match("^(%w+)=")
			if cur then
				state = "properties"
				factories[cur] = {}
			end
		elseif state == "properties" then
			local ps = line:match("properties%((.+)%)")
			if ps then
				props = {}
				for p in ps:gmatch("([%w%.%-_]+)") do
					table.insert(props, p)
				end
				state = "factory"
			end
		elseif state == "factory" then
			if line:find("^/%s*$") then
				state = "new"
			else
				local fac = {}
				table.insert(factories[cur], fac)
				local i = 0
				for p in line:gmatch("([%w%.%-]+)") do
					if i == 0 then
						fac.name = p
					else
						fac[props[i]] = tonumber(p)
					end
					i = i+1
				end
			end
		end
	end
	fp:close()
	return factories
end

local function read_distance(fname, symmetric)
	local mat = {}
	local rows, info = fscg.rows(fname)
	if symmetric then
		local i = 1
		for row in rows, info do
			mat[info.cols[i]] = row
			i = i+1
		end
	else
		for row in rows, info do
			mat[row[info.cols[1]]] = row
		end
	end
	return mat
end

local function parse_expr(expr)
	local terms = {}
	local i = 1
	while i <= #expr do
		local sign = expr:sub(i,i)
		if sign == "+" or sign == "-" then
			i=i+1
		else
			sign = "+"
		end
		local coef
		if string.match(expr:sub(i,i), "[%d%.]") then
			local j = i+1
			while j <= #expr and string.match(expr:sub(j,j), "[%d%.]") do
				j = j+1
			end
			coef = assert(tonumber(expr:sub(i,j-1)))
			if sign == "-" then coef = -coef end
			if expr:sub(j,j) == "*" then j=j+1 end
			i=j
		elseif sign == "-" then
			coef = -1
		end
		local start = i
		while i <= #expr and string.match(expr:sub(i,i), "[^%+%-]") do
			i=i+1
		end
		local var
		if i > start then
			var = expr:sub(start,i-1)
		end
		assert(var or coef)
		table.insert(terms, {var=var, coef=coef})
	end
	return terms
end

local function parse_constraint(src)
	src = string.gsub(src, "%s", "")
	local lhs, operator, rhs = string.match(src, "^([^<>=].+)([<>=])(.+)$")
	local lhs = parse_expr(lhs)
	local rhs = parse_expr(rhs)
	assert(#rhs == 1)
	return {
		lhs = lhs,
		rhs = rhs[1],
		operator = operator
	}
end

local function read_expr(fname)
	local fp = openread(fname)
	local expr = parse_expr(fp:read("*a"))
	fp:close()
	return expr
end

local function read_constraints(fname)
	local fp = openread(fname)
	local buf = buffer.new()
	local constraints = {}
	while true do
		local line = fp:read("*l")
		if not line then break end
		buf:put(line)
		if line:match("[<>=]") then
			table.insert(constraints, parse_constraint(buf:get()))
		end
	end
	return constraints
end

return {
	read_factories   = read_factories,
	read_distance    = read_distance,
	read_expr        = read_expr,
	read_constraints = read_constraints
}
