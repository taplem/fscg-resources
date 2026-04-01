local cda_in = ...

local kuvionumerot = {}
for line in io.open(cda_in, "r"):lines() do
	kuvionumerot[line:match("%d+$")] = true
end

io.stdout:write(io.stdin:read("*L"))

for line in io.stdin:lines("*L") do
	if kuvionumerot[line:match("^%d+")] then
		io.stdout:write(line)
	end
end
