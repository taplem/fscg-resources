local FACTORIES_FILE = "factories.lua"
local UNITS_FILE = "north-karelia.sqlite3"
local OUT_KA_KP_FILE = "pk1_tiem_ka_kp.txt"
local OUT_KP_KP_FILE = "tiem_kp_kp.txt"
local OUT_FACTORIES = "factories_prices.txt"

local sqlite = require "sqlite"
local ffi = require "ffi"

ffi.cdef [[
	typedef struct PJ_AREA PJ_AREA;
	typedef struct PJ_CONTEXT PJ_CONTEXT;
	typedef struct PJ PJ;
	typedef struct { double x, y; } PJ_XY;
	typedef enum proj_direction {
		PJ_FWD   =  1,   /* Forward    */
		PJ_IDENT =  0,   /* Do nothing */
		PJ_INV   = -1    /* Inverse    */
	} PJ_DIRECTION;
	typedef union {
		double v[4];
		PJ_XY xy;
	} PJ_COORD;
	PJ_CONTEXT *proj_context_create(void);
	PJ *proj_create_crs_to_crs(PJ_CONTEXT *ctx, const char *source_crs, const char *target_crs, PJ_AREA *area);
	PJ *proj_normalize_for_visualization(PJ_CONTEXT *ctx, const PJ *obj);
	int proj_trans_array(PJ *P, PJ_DIRECTION direction, size_t n, PJ_COORD *coord);
]]

local proj = ffi.load "proj"
local ctx = proj.proj_context_create()

local tmp_coord = ffi.new("PJ_COORD")
local function transformation(from, to)
	local P = proj.proj_create_crs_to_crs(ctx, from, to, nil)
	return function(x, y)
		tmp_coord.xy.x = x
		tmp_coord.xy.y = y
		assert(proj.proj_trans_array(P, 1, 1, tmp_coord) == 0)
		return tmp_coord.xy.x, tmp_coord.xy.y
	end
end

local wgs84 = transformation("EPSG:2393", "EPSG:4326")
local ykj = transformation("EPSG:4326", "EPSG:2393")

local DEG2RAD = math.pi / 180

print(string.format("[*] reading factories [%s]", FACTORIES_FILE))
local ff = dofile(FACTORIES_FILE)
local ftypes = { "saws", "plywood", "mills", "chips" }
local factories = {}
for _,ft in ipairs(ftypes) do
	for _,f in ipairs(ff[ft]) do
		if f.Y then
			f.type = ft
			f.N_ykj, f.E_ykj = ykj(f.Y, f.X)
			f.E_ykj = f.E_ykj + 3e6
			f.Y = DEG2RAD*f.Y
			f.X = DEG2RAD*f.X
			table.insert(factories, f)
		else
			print(string.format("[!] SKIP %s (missing coordinates)", f.name))
		end
	end
end

local asin, cos, sqrt = math.asin, math.cos, math.sqrt
local function dist(Y1, X1, Y2, X2)
	local hav = 0.5 * (1 - cos(Y2-Y1) + cos(Y1)*cos(Y2)*(1 - cos(X2 - X1)))
	return 1.3 * 6378.137 * 2 * asin(sqrt(hav))
end

---- ka/kp ---------------------------------------------------------------------

local fp_kakp = io.open(OUT_KA_KP_FILE, "w")
fp_kakp:write("standid")
for _,f in ipairs(factories) do
	fp_kakp:write(",", f.name)
end
fp_kakp:write("\n")

print(string.format("[*] reading units [%s]", UNITS_FILE))
local db = sqlite.open(UNITS_FILE)
for row in db:rows("SELECT id, N_ykj, E_ykj FROM stands ORDER BY id") do
	local u = row:row(true)
	local Y, X = wgs84(u.N_ykj*1000, u.E_ykj*1000 + 3e6)
	fp_kakp:write(u.id)
	for _,f in ipairs(factories) do
		fp_kakp:write(" ", dist(f.Y, f.X, DEG2RAD*Y, DEG2RAD*X))
	end
	fp_kakp:write("\n")
end
fp_kakp:close()

---- kp/kp ---------------------------------------------------------------------

local fp_kpkp = io.open(OUT_KP_KP_FILE, "w")
for i,f in ipairs(factories) do
	if i>1 then fp_kpkp:write(",") end
	fp_kpkp:write(f.name)
end
fp_kpkp:write("\n")
for _,f1 in ipairs(factories) do
	for i,f2 in ipairs(factories) do
		if i>1 then fp_kpkp:write(" ") end
		fp_kpkp:write(dist(f1.Y, f1.X, f2.Y, f2.X))
	end
	fp_kpkp:write("\n")
end
fp_kpkp:close()

---- factories_prices ----------------------------------------------------------

-- €/m^3
local PRICE_SS = 80 -- sawlog (softwood)
local PRICE_SH = 70 -- sawlog (hardwood)
local PRICE_PS = 47 -- pulp (softwood)
local PRICE_PH = 49 -- pulp (hardwood)
local PRICE_EW = 36 -- energy wood

local sawmills = {}
local pulpmills = {}
local chips = {}

for _,f in ipairs(factories) do
	if f.type == "saws" or f.type == "plywood" then
		table.insert(sawmills, f)
	elseif f.type == "mills" then
		table.insert(pulpmills, f)
	else
		table.insert(chips, f)
	end
end

local fp_fac = io.open(OUT_FACTORIES, "w")
fp_fac:write("sawmill=list(")
for i,f in ipairs(sawmills) do
	if i>1 then fp_fac:write(",") end
	fp_fac:write(f.name)
end
fp_fac:write(")\n")
fp_fac:write("properties(xcoor,ycoor,scap1,sp1,sprice1,scap2,sp2,sprice2,scapsoftw,scaptot,scap3,sp3,sprice3,bark1,bark2,bark3,dust,chip)\n")
for _,f in ipairs(sawmills) do
	fp_fac:write(
		f.name,
		",", f.E_ykj, -- xcoor
		",", f.N_ykj, -- ycoor
		",", (f.pine or 0)*f.technical_cap, -- scap1
		",1", -- sp1
		",", PRICE_SS, -- sprice1
		",", (f.spruce or 0)*f.technical_cap, -- scap2
		",2", -- sp2
		",", PRICE_SS, -- sprice2
		",", math.min((f.pine or 0)+(f.spruce or 0), 1)*f.technical_cap, -- scapsoftw
		",", f.technical_cap, -- scaptot
		",", (f.birch or 0)*f.technical_cap, -- scap3
		",3", -- sp3
		",", PRICE_SH, -- sprice3
		",0", -- bark1 (???)
		",0", -- bark2 (???)
		",0", -- bark3 (???)
		",0", -- dust (???)
		",", f.type == "saws" and 0.3 or 0.315, -- chips
		"\n"
	)
end
fp_fac:write("/\n")

fp_fac:write("pulpmill=list(")
for i,f in ipairs(pulpmills) do
	if i>1 then fp_fac:write(",") end
	fp_fac:write(f.name)
end
fp_fac:write(")\n")
fp_fac:write("properties(xcoor,ycoor,pcap1,sp1,pprice1,pcap2,sp2,pprice2,pcapsoftw,pcaptot,pcap3,sp3,pprice3,mixsoftw)\n")
for _,f in ipairs(pulpmills) do
	fp_fac:write(
		f.name,
		",", f.E_ykj, -- xcoor
		",", f.N_ykj, -- ycoor
		",", (f.pine or 0)*f.technical_cap, -- pcap1
		",1", -- sp1
		",", PRICE_PS, -- pprice1
		",", (f.spruce or 0)*f.technical_cap, -- pcap2
		",2", -- sp2
		",", PRICE_PS, -- pprice2
		",", math.min((f.pine or 0)+(f.spruce or 0), 1)*f.technical_cap, -- pcapsoftw
		",", f.technical_cap, -- pcaptot
		",", (f.birch or 0)*f.technical_cap, -- pcap3
		",3", -- sp3
		",", PRICE_PH, -- pprice3
		",0", -- mixsoftw (???)
		"\n"
	)
end
fp_fac:write("/\n")

fp_fac:write("chips=list(")
for i,f in ipairs(chips) do
	if i>1 then fp_fac:write(",") end
	fp_fac:write(f.name)
end
fp_fac:write(")\n")
fp_fac:write("properties(xcoor,ycoor,meha_im3,meha_kim3,teolha_kim3,hakeYht_kim3,cprice,cprice2,karsittu_p,kokopuu_p,hakkuutahde,jarea_p,kannot_p,kuitu_p,muu_p)\n")
for _,f in ipairs(chips) do
	fp_fac:write(
		f.name,
		",", f.E_ykj, -- xcoor
		",", f.N_ykj, -- ycoor
		",0", -- meha_im3 (???)
		",", f.technical_cap, -- meha_kim3
		",", f.technical_cap, -- teolha_kim3
		",0", -- hakeYht_kim3 (???)
		",", PRICE_EW, -- cprice
		",", PRICE_EW, -- cprice2
		",0", -- karsittu_p (???)
		",0", -- kokopuu_p (???)
		",0", -- hakkuutahde (???)
		",0", -- jarea_p (???)
		",0", -- kannot_p (???)
		",0", -- kuitu_p (???)
		",0", -- muu_p (???)
		"\n"
	)
end
fp_fac:write("/\n")

fp_fac:write [[
homeuse=list(ppESKT)
properties(xcoor,ycoor,captot)
ppESKT,6840892,512009,315000
/
]]

fp_fac:close()
