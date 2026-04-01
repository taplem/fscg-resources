require "vmi"
local lim3 = require "lim3"

local function site(x)
	return string.format("site._'{%s}", x)
end

-- user must define site.tending_cf
-- (TODO: use query param, see comment in select_ba_thinning)
local select_tending = lim3.selector {
	var = "f",
	select = "F - tending_cf*mh_regeneration_f1", -- tpoist = ttot - rr*tp(3)
	{
		which = "all(mc <= 1, d <= 18)",
		priority = "v",
		-- y = {0.5, 0},
		y = {50, 0},
		method = "level"
	}
}

local function tending_selector(cf)
	return lim3.wrap(select_tending):config("site", {tending_cf=cf})
end

local select_first_thinning = lim3.selector {
	var = "f",
	select = "F - mh_first_thinning_f1",
	{
		which = "all(mc <= 1, d >= 15)",
		priority = "d",
		-- y = {0, 1},
		y = {0, 100},
		method = "level"
	},
	{
		which = "all(mc <= 1, s != site.sdom)",
		priority = "v",
		-- y = {0.5, 0},
		y = {50, 0},
		method = "level"
	},
	{
		which = "mc <= 1",
		priority = "v",
		-- y = {0.5, 0},
		y = {50, 0},
		method = "level"
	}
}

-- user must define site.ba_thinning_dmin
-- (TODO: use query.ba_thinning_dmin here instead, need query param support in lim3.config)
local select_ba_thinning = lim3.selector {
	var = "g",
	-- TODO: target should be min(0.5*G, G-mh_thinning_g1) ???
	select = "G - mh_thinning_g1",
	{
		which = "all(d >= site.ba_thinning_dmin, mc <= 1, s != site.sdom)",
		priority = "v",
		-- y = {0.5, 0},
		y = {50, 0},
		method = "level"
	},
	{
		which = "all(d >= site.ba_thinning_dmin, mc <= 1)",
		priority = "v",
		-- y = {0.5, 0},
		y = {50, 0},
		method = "level"
	}
}

local function ba_thinning_selector(dmin)
	return lim3.wrap(select_ba_thinning):config("site", {ba_thinning_dmin=dmin})
end

local select_clearcut = lim3.selector {
	which = "all(d >= 6, mc <= 1)"
}

local select_seedtree_cutting = lim3.selector {
	var = "f",
	select = "F - mh_retention_f1",
	{
		which = "all(d >= 6, mc <= 1, s != 1)",
	},
	{
		which = "all(d >= 6, mc <= 1)",
		priority = "h",
		y = {5, 0},
		method = "level"
	}
}

local check_energy = control.check "site.h_RE >= 15*site.area"
local check_energy0 = control.check "site.h_RE > 0"
local check_first_thinning_yield = control.check "site.h_F >= 500*site.area"

local thin_e0 = lim3.cutting()
local thin_e4 = lim3.cutting(check_energy):config("site", { p_bewtag=4 })
local thin_e6 = lim3.cutting(check_energy):config("site", { p_bewtag=6 })
local thin_e7 = lim3.cutting(check_energy):config("site", { p_bewtag=7 })
local thin_e9 = lim3.cutting(check_energy):config("site", { p_bewtag=9 })

local cc_e0   = lim3.cutting():config("site", { p_btag=2 })
local cc_e1   = lim3.cutting(check_energy0):config("site", { p_btag=2, p_bewtag=1 })
local cc_e3   = lim3.cutting(check_energy0):config("site", { p_btag=2, p_bewtag=3 })

local function flag(event)
	return data.transaction()
		:update("site", { [string.format("flag'{%s}", event)] = "year" })
end

local exec_retention = control.all {
	lim3.selector {
		var = "v",
		select = 5,
		{
			which = "all(d >= 7, mc <= 1)",
			priority = "v",
			y = {0, 0.5},
			method = "level"
		},
	},
	lim3.selector {
		var = "w",
		select = 100
	},
	lim3.split { mc=2 }
}

local exec_clear_regeneration_area = control.all {
	lim3.selector {
		{ which = "all(d<=15, any(s<=2, s=5, s=9))" },
		-- TODO: { which="any(d>15, s=3, s=4, s=6, s=7, s=8)", retain=1500 }
		-- (see comments in select.lua)
	},
	-- TODO: this isn't really a thinning for cost purposes.
	--       cutting costs should be counted as regeneration costs.
	-- (the most efficient implementation is probably a version of getcut() for each different
	--  type of cutting for which costs need to be counted separately + one shared getcut() that
	--  computes all the common variables)
	thin_e0
}

local try_site_preparation = control.all {
	control.try(control.all {
		site "all(mty<=5, mc<=2)",
		exec_clear_regeneration_area,
		lim3.trace "35 CLEAR REGEN"
	}),
	control.try(control.all {
		site "all(mty>=2, mty<=4, mc<=1)",
		lim3.soilprep(),
		lim3.trace "36 SITE PREP"
	})
}

local events = control.all {
	-- First thinnings ---------------------
	control.try(control.all {
		site [[
			max(
				flag'first_thinning_blocked,
				flag'basal_area_thinning,
				flag'tending
			) <= year-10
			and F'[mc<=1] >= 1500
			and dg'[mc<=1] >= 8.0
			and hg'[mc<=1] >= 5.0
			and hg'[mc<=1] <= 13.5
		]],
		control.optional(control.all {
			select_first_thinning,
			check_first_thinning_yield,
			control.any {
				-- 10 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				control.all {
					site [[ any(all(mc<=2, alr=1), all(mc<=1, any(alr=2, alr=3), t_drain<=1950),
						all(mc=2, alr>=2)) ]],
					thin_e0,
					lim3.trace "10 FIRST THINNING",
				},
				-- 11 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				-- TODO: min_interval: 20
				control.all {
					site "all(mc<=1, t_drain>1950, any(alr>2, all(alr=2, mty>=4, mty<=6)))",
					thin_e0,
					lim3.trace "11 FIRST THINNING",
					-- TODO: DRAINING type=2
				},
				-- kommentoitu pois MS_EVENTissä
				-- -- 12 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				-- control.all {
				-- 	site "all(mc<=1, alr=1, mty<=4, sdom!=2)",
				-- 	thin_e4,
				-- 	lim3.trace "12 FIRST THINNING",
				-- },
				-- 13 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				control.all {
					site [[ all(mc<=1, da'[mc<=1]<20,
						any(alr=1, all(t_drain<=1950, any(alr=2, alr=3)),
							all(t_drain>1950, alr=2, mty<=3)))]],
					thin_e6,
					lim3.trace "13 FIRST THINNING",
				},
				-- 14 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				-- TODO: min_interval: 20
				control.all {
					site "all(mc<=1, da'[mc<=1]<20, t_drain>1950, any(all(alr=2, mty>=4), alr>=3))",
					thin_e6,
					-- TODO: DRAINING type=2
					lim3.trace "14 FIRST THINNING",
				},
				-- 15 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				control.all {
					site "all(mc<=1, alr=1, mty<=4, sdom!=2, da'[mc<=1]<15) ",
					thin_e7,
					lim3.trace "15 FIRST THINNING",
				},
				-- 16 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				control.all {
					site "all(mc<=1, da'[mc<=1]<18, any(alr=1, all(alr=2, mty<=3, t_drain>1950)))",
					thin_e9,
					lim3.trace "16 FIRST THINNING",
				},
				-- 17 FIRST THINNING - NUMBER OF STEMS/HA INSTRUCTIONS
				control.all {
					site [[ all(mc<=1, t_drain>1950, da'[mc<=1]<18,
						any(all(alr=2, mty>=4, mty<=6), alr>=3)) ]],
					thin_e9,
					-- TODO: DRAINING type=2
					lim3.trace "17 FIRST THINNING",
				}
			},
			flag "first_thinning"
		}),
		flag "first_thinning_blocked",
	}),

	-- Basal area thinnings ----------------
	control.try(control.all {
		site [[
			max(
				flag'basal_area_thinning_blocked,
				flag'first_thinning,
				flag'seedtree_cutting,
				flag'clearcut
			) <= year-10
			and G'[mc<=1] >= 0.9*mh_thinning_g0
			and hg'[mc<=1] >= 13.5
			and dg'[mc<=1] <= 1.5*mh_regeneration_d0
			and ag'[mc<=1] <= 1.5*mh_regeneration_a0
		]],
		control.optional(control.all {
			control.any {
				control.all {
					ba_thinning_selector(6),
					control.any {
						-- 18 THINNING - BASAL AREA INSTRUCTIONS
						control.all {
							site [[ any(all(mc<=2, alr=1), all(mc<=1, any(alr=2, alr=3), t_drain<=1950),
									all(mc<=1, alr=2, mty<=3, t_drain>1950), all(mc=2, alr>=2)) ]],
							thin_e0,
							lim3.trace "18 THINNING",
						},
						-- 19 THINNING - BASAL AREA INSTRUCTIONS + OJITUS
						-- TODO: min_interval: 20
						control.all {
							site [[ all(mc<=1, t_drain>1950,
								any(all(alr=2, mty>=4, mty<=6), all(alr>=3, alr<=4))) ]],
							thin_e0,
							-- TODO: DRAINING type=2
							lim3.trace "19 THINNING",
						}
					}
				},
				control.all {
					site "mc<=1",
					ba_thinning_selector(4),
					control.any {
						-- kommentoitu pois MS_EVENTissä
						-- -- 20 THINNING - BASAL AREA INSTRUCTIONS
						-- control.all {
						-- 	site "all(alr=1, mty<=4, sdom!=2)",
						-- 	thin_e4,
						-- 	lim3.trace "20 THINNING",
						-- },
						-- 21 THINNING - BASAL AREA INSTRUCTIONS
						control.all {
							site "all(da'[mc<=1]<20, any(alr=1, all(alr=2, mty<=3, t_drain>1950)))",
							thin_e6,
							lim3.trace "21 THINNING",
						},
						-- 22 THINNING - BASAL AREA INSTRUCTIONS
						-- TODO: min_interval: 20
						control.all {
							site [[ all(t_drain>1950, da'[mc<=1]<20,
								any(all(alr=2, mty>=4, mty<=6), all(alr>=3, alr<=4))) ]],
							thin_e6,
							-- TODO: DRAINING type=2
							lim3.trace "22 THINNING",
						},
						-- 23 THINNING - BASAL AREA INSTRUCTIONS
						control.all {
							site "all(alr=1, mty<=4, sdom!=2, da'[mc<=1]<15) ",
							thin_e7,
							lim3.trace "23 THINNING",
						},
						-- 24 THINNING - BASAL AREA INSTRUCTIONS
						control.all {
							site "all(da'[mc<=1]<18, any(alr=1, all(alr=2, mty<=3, t_drain>1950)))",
							thin_e9,
							lim3.trace "24 THINNING",
						},
						-- 25 THINNING - BASAL AREA INSTRUCTIONS
						-- TODO: min_interval: 20
						control.all {
							site [[ all(t_drain>1950, da'[mc<=1]<18,
								any(all(alr=2, mty>=4, mty<=6), all(alr>=3, alr<=3))) ]],
							thin_e9,
							-- TODO: DRAINING type=2
							lim3.trace "25 THINNING",
						}
					}
				}
			},
			flag "basal_area_thinning"
		}),
		flag "basal_area_thinning_blocked"
	}),

	-- Over story removal (TODO) -----------
	-- 26 OVER STORY REMOVAL
	-- * mc in {0,1,2}
	-- comparable: 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 30, 31, 32
	-- precedessors: 99
	-- * CUTTING(4.0) grem=1.0, grminu=0.5, grminp=2.0, drmin=0.7, amin=0.85, frmin=0.8, season=2.0, logm=0.0, logsk=0.0

	-- Seed tree and shelterwood cutting ---
	control.try(control.all {
		site [[
			max(
				flag'seedtree_cutting_blocked,
				flag'first_thinning,
				flag'basal_area_thinning,
				flag'clearcut
			) <= year-10
			and G'[mc<=1] >= 0.7*mh_thinning_g1
			and ag'[mc<=1] <= 1.6*mh_regeneration_a0
			and (
				dg'[mc<=1] >= 0.95*mh_regeneration_d0
				or ag'[mc<=1] >= 0.95*mh_regeneration_a0
			)
		]],
		control.optional(control.all {
			control.any {
				-- 27 SEED TREE CUTTING FOR NATURAL REGENERATION (PINE)
				control.all {
					site "all(mc<=1, sdom=1, any(all(mty>=4, mty<=6), all(mty>=1, mty<=3, dd<=750)))",
					exec_retention,
					select_seedtree_cutting,
					thin_e0,
					lim3.trace "27 SEED TREE CUTTING",
				},
				-- 29 SHELTERWOOD CUTTING FOR NATURAL REGENERATION (SPRUCE)
				control.all {
					site "all(mc=2, mty<=3, sdom=2)",
					exec_retention,
					select_clearcut,
					thin_e0,
					lim3.trace "29 SHELTERWOOD CUTTING",
				}
			},
			flag "seedtree_cutting",
			try_site_preparation
		}),
		flag "seedtree_cutting_blocked"
	}),

	-- Clear cutting -----------------------
	control.try(control.all {
		site [[
			max(
				flag'clearcut_blocked,
				flag'first_thinning,
				flag'basal_area_thinning,
				flag'seedtree_cutting
			) <= year-10
			and (
				dg'[mc<=1] >= 0.95*mh_regeneration_d0
				or ag'[mc<=1] >= select(mty<6, 0.95, 2.0)*mh_regeneration_a0
				or (dg'[mc<=1] >= 8 and G'[mc<=1] < 0.5*mh_thinning_g1)
				or (dg'[mc<=1] < 8 and F'[mc<=1] < mh_young_stand_fmin)
			)
		]],
		control.optional(control.all {
			exec_retention,
			select_clearcut,
			control.any {
				-- 30 CLEAR CUTTING
				control.all {
					site "all(mc<=1, mty<=4)",
					cc_e0,
					lim3.trace "30 CLEAR CUTTING",
				},
				-- 31 CLEAR CUTTING
				control.all {
					site "all(mc<=1, alr=1, mty>=2, mty<=4, verlt=0)",
					cc_e1,
					lim3.trace "31 CLEAR CUTTING",
				},
				-- 32 CLEAR CUTTING
				control.all {
					site "all(mc<=1, alr=1, mty>=2, mty<=3, verlt=0, sdom=2)",
					cc_e3,
					lim3.trace "32 CLEAR CUTTING",
				}
			},
			flag "clearcut",
			try_site_preparation
		}),
		flag "clearcut_blocked"
	}),

	-- Regeneration ------------------------
	control.try(control.all {
		site [[
			max(
				flag'regeneration_blocked,
				flag'first_thinning,
				flag'basal_area_thinning,
				flag'seedtree_cutting,
			) <= year-10
			and mc <= 1
			# and year <= 2065 # (testi)
			and G'[mc<=1] < 8
			and (
				F'[mc<=1] <= 10
				or all(flag'regeneration<flag'clearcut, year-flag'clearcut<10)
				or all(
					mty <= 3,
					isdecid(sdom),
					any(F'[mc<=1] < 900, year-flag'regeneration > 30),
					dg'[mc<=1] >= 2,
					G'[mc<=1] > 1
				)
			)
		]],
		control.any {
			control.all {
				control.any {
					-- 41 REGENERATION OF PINE
					control.all {
						site "all(mty>=3, mty<=4)",
						-- TODO: uusi viljely: ota puun tiedot ja generoi hajontaa.
						lim3.regeneration {
							s = 1,
							-- TODO: fhk linspace(), tai vielä parempi, (vektoriarvoinen) rand()
							h = {0.050, 0.075, 0.100, 0.125, 0.150},
							a = {0.50,  0.75,  1.00,  1.25,  1.50},
							snt = 2,
						},
						lim3.trace "41 REGENERATION PINE",
					},
					-- 42 REGENERATION OF SPRUCE
					control.all {
						site "all(mty>=2, mty<=3)",
						lim3.regeneration {
							s = 2,
							h = {0.10, 0.15, 0.20, 0.25, 0.30},
							a = {1.50, 1.75, 2.00, 2.25, 2.50},
							snt = 2,
						},
						lim3.trace "42 REGENERATION SPRUCE",
					},
					-- 43 REGENERATION OF BIRCH
					control.all {
						site "all(mty>=1, mty<=3)",
						lim3.regeneration {
							s = 3,
							h = {0.30, 0.35, 0.40, 0.45, 0.50},
							a = {0.50, 0.75, 1.00, 1.25, 1.50},
							snt = 2,
						},
						lim3.trace "43 REGENERATION BIRCH",
					}
				},
				flag "regeneration"
			},
			-- ei sallita ohitusta jos tehtiin avohakkuu
			-- (TODO: pakotettu täydennysistutus korvaa tämän ehdon)
			control.check "site.flag'clearcut < site.year"
		},
		flag "regeneration_blocked",
	}),

	-- Supplementary planting --------------
	-- TODO: poissa käytöstä nyt koska satunnaistus puuttuu
	-- 45 SUPPLEMENTARY PLANTING
	-- branching: 0
	-- * mc in {0,1,2}
	-- precedessors: 26, 41, 42, 43, 99
	-- * REGENERATION type=3.0, years=3.0, spe=0.0, f=0.0, a=2.0, h=0.3, sp=0.0, mtymin=1.0,
	--   mtymax=5.0
	-- local event45 = control.all {
	-- 	"mc <= 2 and mty >= 1 and mty <= 5",
	-- 	lim3.planting {
	-- 
	-- 	}
	-- }
	-- (TODO: satunnaistus puuttuu)
	-- control.try(control.all {
	-- 	check [[
	-- 		any(
	-- 			flag'regeneration > year-10,
	-- 			max(
	-- 				flag'first_thinning,
	-- 				flag'basal_area_thinning,
	-- 				flag'clearcut
	-- 			) <= year-10
	-- 		)
	-- 	]],
	-- 	event45
	-- })

	-- Pre-tending -------------------------
	control.try(control.all {
		site [[
			max(
				flag'first_thinning,
				flag'pre_tending,
				flag'tending
			) <= year-10
			and mc <= 2
			and F'[mc<=1] >= max(1000, 2.5*mh_regeneration_f1)
			and all(dg'[mc<=1] >= 0.5, dg'[mc<=1] <= 3.5)
			and all(hg'[mc<=1] >= 0.75, hg'[mc<=1] <= 2.0)
		]],
		tending_selector(1.5),
		thin_e0,
		lim3.trace "49 PRE-TENDING",
		flag "pre_tending"
	}),

	-- Tending of young stands -------------
	-- TODO: koodissa hardkoodattuna ehto, ettei tätä tehdä jos edellinen hakkuu oli harvennushakkuu,
	--       tässä pitäisi varmaan myös olla vastaava ehto
	control.try(control.all {
		site [[
			max(
				flag'first_thinning,
				flag'basal_area_thinning,
				flag'tending
			) <= year-10
			and mc <= 2
			and F'[mc<=1] >= max(200, 1.15*mh_regeneration_f1)
			and all(dg'[mc<=1] >= 1.0, dg'[mc<=1] <= 7.7)
			and all(hg'[mc<=1] >= 2.0, hg'[mc<=1] <= 6.0)
		]],
		tending_selector(1.0),
		thin_e0,
		lim3.trace "50 TENDING",
		flag "tending"
	})
}

--------------------------------------------------------------------------------

local SAMPLE_RATIO = 1       -- joka n's kuvio
local DATA_COEF = 1.1      -- metsäkeskusaineisto -> koko pohjois-karjala
local ACOEF = DATA_COEF*SAMPLE_RATIO  -- otos -> koko pohjois-karjala
local PCOEF = ACOEF/10                -- kausittainen (10v) -> vuosittainen

data.define [[
	model site value = area * sum(tree.f*tree.value_t)
	model site h_NUt4 = select(year<=2075, 0, h_NU*1.04^(2075-year))
]]

local nodes = {
	2035, 2045, 2055, 2065, 2075,
	["report.nodes"] = {
		stand = "site.id",
		area = "site.area",
		year = "site.year",
		RC = PCOEF.."*site.RC",
		RC1 = PCOEF.."*site.RC'[s=1]",
		RC2 = PCOEF.."*site.RC'[s=2]",
		RCLP = PCOEF.."*site.RC'[any(all(s>=5,s<=7),s=9)]",
		RL = PCOEF.."*site.RL",
		RL1 = PCOEF.."*site.RL'[s=1]",
		RL2 = PCOEF.."*site.RL'[s=2]",
		RL3 = PCOEF.."*site.RL'[any(s=3,s=4)]",
		RP = PCOEF.."*site.RP",
		RP1 = PCOEF.."*site.RP'[s=1]",
		RP2 = PCOEF.."*site.RP'[s=2]",
		RP3 = PCOEF.."*site.RP'[any(s=3,s=4)]",
		RE = PCOEF.."*site.RE",
		REst = PCOEF.."*site.REst",
		REcr = PCOEF.."*site.REcr",
		REsr = PCOEF.."*site.REsr",
		WC = PCOEF.."*site.WC",
		WEc = PCOEF.."*site.WEc",
		NU = PCOEF.."*site.NU",
		WRA = PCOEF.."*site.WRA",
	},
}

local tail = {
	-- site "max(flag'clearcut, flag'seedtree_cutting) = year",
	site "max(flag'clearcut, flag'seedtree_cutting) = year or year >= 2200",
	["report.tail"] = {
		stand = "site.id",
		year = "site.year",
		NUt4 = ACOEF.."*site.NUt4",
		-- value = AREA_COEF.."*site.value",
		ag = "site.ag",
		hg = "site.hg",
		dg = "site.dg",
		G = "site.G",
		F = "site.F"
	}
}

--------------------------------------------------------------------------------

local dbin, dbout = ...
data.attach(dbin or error("no input data specified"))
if dbout then
	data.attach(dbout, "report")
end

lim3.setup {
	output = dbout and "db" or "console",
	grow = "metsi",
	events = events,
	nodes = nodes,
	tail = tail,
}

data.task = string.format("SELECT id FROM stands WHERE id%%%d=0", SAMPLE_RATIO)
