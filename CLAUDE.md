# CLAUDE.md — Astro Magick

Developer notes for `astro_magick.py`.

---

## Overview

Command-line astronomical almanac for magickal practice.
Powered by **Skyfield + JPL DE421 ephemeris** (~17 MB, downloaded on first run).
Behenian fixed star data requires the **Hipparcos catalog** (`hip_main.dat`, ~37 MB,
downloaded on first `--behenian` run).

```
python3 astro_magick.py
python3 astro_magick.py --date 2026-06-21
python3 astro_magick.py --lat 51.5074 --lon -0.1278 --timezone Europe/London
python3 astro_magick.py --days 7
python3 astro_magick.py --nighttransit
python3 astro_magick.py --behenian
python3 astro_magick.py --behenian --nighttransit
python3 astro_magick.py --eclipses
python3 astro_magick.py --json | jq .
```

---

## Dependencies

```
pip install skyfield
```

DE421 covers 1900–2053-10-09. Eclipse scan hard-stops at 2053-10-01.
For dates outside that range, switch to `de440.bsp` in `load_ephemeris()`.

---

## Configuration file

`astro_magick.cfg` is read from the current directory. Created with defaults
on first run if absent.

```ini
[location]
latitude  = 37.3382
longitude = -121.8863
name      = San Jose, CA

[time]
timezone = America/Los_Angeles
```

CLI flags `--lat`, `--lon`, `--timezone` override the config file.

---

## Command-line flags

  --lat           Latitude (decimal degrees, + = North)
  --lon           Longitude (decimal degrees, + = East)
  --timezone      IANA timezone name
  --date          Start date YYYY-MM-DD (default: today)
  --days N        Number of consecutive daily reports (default: 1)
  --nighttransit  Next night-time peak for each planet + Behenian star peaks (if --behenian)
  --behenian      Behenian fixed star aspects; with --nighttransit adds star peak table
  --eclipses      Next partial and total lunar/solar eclipses (up to 2053)
  --json          All requested reports as a single JSON object to stdout

All status/error messages go to stderr; stdout is clean for piping.
There is no --astronomical flag; IAU constellation data appears in the main table.

---

## POSITIONS table

Printed once per day (respects --days). Columns:

  BODY           Planet name
  SIGN           Tropical zodiac sign
  DEG            Degrees within sign (0.0–29.9)
  Rx             " Rx" if retrograde, else blank
  CONSTELLATION  IAU constellation + % through it (e.g. "Gemini 72%")
  RISE           Local time of horizon rise
  TRANSIT        Upper culmination. [D] = Sun up, [N] = Sun down
  SET            Local time of horizon set
  ANTITRANSIT    Lower culmination (nadir). [D]/[N] as above
  NVIS           Start of night visibility window
  END            End of night visibility window

Notes:
- CONSTELLATION column sits between Rx and RISE
- Longest constellation string: "Sagittarius 100%" = 16 chars; column width = 17
- [D]/[N] tags are always 3 ASCII chars — no emoji
- Circumpolar bodies show `----` for RISE/SET
- All column content is pure ASCII; planet glyphs are non-padded decorative prefixes

### IAU ecliptic boundary data

Source: David Asher, Human Orrery project (cantab.net), J2000 ecliptic longitudes.
The ecliptic crosses 13 constellations (Ophiuchus at 247.6–266.2°):

  Pisces        0.000      Aries        28.687     Taurus       53.417
  Gemini       90.140      Cancer      117.988     Leo         138.038
  Virgo       173.851      Libra       217.810     Scorpius    241.047
  Ophiuchus   247.638      Sagittarius 266.238     Capricornus 299.656
  Aquarius    327.488      Pisces      351.650 (resumes)

---

## NEXT NIGHT PEAK table  (--nighttransit)

Runs once per day in the --days loop (each day scans forward from that date).
Columns: BODY | DATE | TIME | SIGN | DEG | CONSTELLATION | TYPE

### Moon rows (3 per run)

  Moon       Next night any phase is above the horizon (elongation 0–360°)
  Full Moon  Next Full Moon night (elongation 165–195°)
  New Moon   Next New Moon night (elongation 345–15°, wrapping window)

All three use `find_next_night_peak_moon_phase()`.

### Other bodies

  Sun               Always daytime — shown as `---`
  Mercury, Venus    Inner planets: 30-min sampling across night window
  Mars, Jupiter,    `find_next_night_transit_outer()`: enumerates all transits
  Saturn            in one Skyfield call; returns first inside a night window

### TYPE column

  transit     Meridian crossing falls within the night window
  night peak  Transit is in daytime; shown value is highest point during darkness

---

## BEHENIAN FIXED STARS  (--behenian)

### Daily aspects table

Runs once per day in the --days loop. Shows every planet that is within
BEHENIAN_WARN (20°) of a Behenian star AND approaching it (i.e. Skyfield
finds a closer approach in the future). Rows where the planet is moving away
are silently dropped.

  PLANET    Body name
  STAR      Behenian star name
  SEP       Current ecliptic separation in degrees
  STATUS    "IN ORB  <6 deg" (within traditional 6° orb) or "approaching"
  APPROACH  Local datetime of closest approach (Skyfield-computed, not estimated)
  VISIBLE   "night" if star is >15° from Sun; "solar glare" otherwise

### Conjunction time algorithm  (`find_conjunction_time`)

Previous approach (find_discrete + 3° threshold) was wrong: if a planet
approaches to only 6° it never crosses 3°, so find_discrete saw nothing and
returned None → row was incorrectly dropped or labelled "none in 1 yr".

Current approach:
1. Sample ecliptic separation every 6 hours across the full 366-day window
2. Find the global minimum separation across all samples
3. If minimum ≥ current separation (planet moving away) → return None → row dropped
4. Otherwise ternary-refine the minimum to ~second accuracy (60 iterations, ±2 days)
5. Return the local datetime of closest approach

This correctly handles all cases: prograde approach, retrograde approach,
close passes that never reach a fixed threshold (like Mars/Fomalhaut at 5°),
and bodies that are moving away.

### Star catalog  (16 stars)

All stars identified by Hipparcos catalog number, loaded via
`skyfield.data.hipparcos` + `skyfield.api.Star`. No hand-coded positions.
Catalog cached in `_HIP_DF_CACHE` / `_BEHENIAN_STAR_CACHE` (loaded once per run).

  Algol         HIP 14576   Beta Persei
  Alcyone       HIP 17702   Eta Tauri (Pleiades)
  Aldebaran     HIP 21421   Alpha Tauri
  Capella       HIP 24608   Alpha Aurigae
  Sirius        HIP 32349   Alpha Canis Majoris
  Procyon       HIP 37279   Alpha Canis Minoris
  Regulus       HIP 49669   Alpha Leonis
  Alkaid        HIP 67301   Eta Ursae Majoris
  Algorab       HIP 59803   Delta Corvi
  Spica         HIP 65474   Alpha Virginis
  Arcturus      HIP 69673   Alpha Bootis
  Alphecca      HIP 76267   Alpha Coronae Borealis
  Antares       HIP 80763   Alpha Scorpii
  Vega          HIP 91262   Alpha Lyrae
  Deneb Algedi  HIP 107556  Delta Capricorni
  Fomalhaut     HIP 113368  Alpha Piscis Austrini

Multiple rows for one planet are normal and correct — each row is a distinct
planet×star aspect. The Sun near both Deneb Algedi and Fomalhaut (which are
~33° apart) will appear twice.

---

## BEHENIAN STARS NIGHT PEAK table  (--behenian --nighttransit)

Runs once per day in the --days loop. For each of the 16 Behenian stars:

  STAR               Star name
  DATE               Next date star is above horizon at night
  TIME               Local time of highest point during that night
  TYPE               "transit" if meridian crossing falls in night window;
                     "night peak" if transit is in daytime
  NEXT NIGHT TRANSIT Next date+time the star's meridian crossing falls in darkness
                     (scans up to 2 years)
  CONJ               Any planet within 6° of the star at the peak moment

---

## --eclipses table

Scans up to 2053-10-01 (DE421 limit) for:

  Lunar  partial   Next penumbral or partial lunar eclipse
  Lunar  total     Next total lunar eclipse
  Solar  partial   Next partial solar eclipse visible at location
  Solar  total     Next total or annular solar eclipse at location

Lunar: `skyfield.eclipselib.lunar_eclipses()`.
Solar: topocentric Sun–Moon separation scan + ternary refinement.
All ephemeris calls wrapped in try/except for DE421 boundary safety.

---

## --json output

All reports packed into one JSON object. Each day in `"almanac"` contains
the positions, and optionally `"behenian"` aspects and `"night_transit"` data
for that specific day (respects --days correctly).

  out["almanac"]          always; list of one entry per --days day
    entry["behenian"]     --behenian aspects for that day
    entry["night_transit"] --nighttransit peaks scanning from that day
    entry["behenian_night"] --behenian + --nighttransit star peaks from that day
  out["eclipses"]         --eclipses

---

## Alignment

Column widths defined in dicts (W, NW, BW, SNW) shared between header and all
data rows via a single fmt_* closure — misalignment is structurally impossible.
All column content is pure ASCII. Planet/zodiac glyphs are non-padded prefixes.

---

## Key functions

  astronomical_constellation()       IAU boundary lookup → (name, pct_through)
  behenian_aspects()                 Planet→star aspects; drops moving-away rows
  behenian_ecl_lon()                 Ecliptic lon of a Star object at time t
  build_behenian_json()              JSON for --behenian daily aspects
  build_behenian_night_json()        JSON for --behenian --nighttransit star peaks
  build_daily_json()                 JSON for one day's almanac entry
  build_eclipses_json()              JSON for --eclipses
  build_nighttransit_json()          JSON for --nighttransit (one day)
  check_conjunctions()               Planet pairs within threshold degrees
  check_eclipses()                   Same-day eclipse warning for daily header
  ecl_separation()                   Shortest angular distance (0–180°)
  find_conjunction_time()            Sampling-based actual conjunction datetime
  find_eclipses()                    Scan for lunar + solar eclipses
  find_events()                      Rise, upper transit, set, lower transit
  find_next_night_peak()             Night peak dispatcher for planets
  find_next_night_peak_moon_phase()  Moon night peak for arbitrary elongation range
  find_next_night_transit_outer()    Bulk transit scan for Mars/Jupiter/Saturn
  find_next_star_night_peak()        Night peak for a Behenian Star object
  find_next_star_night_transit()     Next nighttime transit for a Star object
  find_star_night_peak()             Highest point within one night window
  find_star_transit()                Meridian transit within a time window
  get_behenian_stars()               Load/cache Hipparcos Star objects (16 stars)
  get_body_position()                Ecliptic lon, RA, Dec at given time
  get_moon_phase()                   Phase name + illumination % from elongation
  is_retrograde()                    Detect retrograde by sampling ecl lon ±6h
  load_config()                      Read astro_magick.cfg, create if missing
  load_ephemeris()                   Load DE421, download if missing
  night_peak_altitude()              30-min sampling + ternary refinement
  night_window()                     Compute sunset → next-sunrise window
  output_json()                      Collect all reports into one JSON object
  planet_night_visibility()          Clip planet rise-set to night window
  print_behenian_night_report()      --behenian --nighttransit star peak table
  print_behenian_report()            --behenian daily aspects table
  print_eclipse_report()             --eclipses table
  print_night_transit_report()       --nighttransit planet table
  print_report()                     Daily positions table
  star_alt_at()                      Altitude of a Star object at given TT JD
  zodiac_sign()                      Tropical sign from ecliptic longitude
  _get_night_window_for_day()        Sunset→sunrise helper for one calendar day
  _utc_in_night()                    Check if UTC datetime falls in night window

Note: `print_astronomical_report`, `print_astronomical_night_report`,
`build_astronomical_json`, and `build_astronomical_night_json` remain in the
file as dead code from the pre-merge era. They are not called anywhere and
can be removed in a future cleanup pass.

---

## Accuracy

  Positions:    ~arcsecond level (DE421 + Skyfield apparent place)
  Fixed stars:  Hipparcos catalog positions propagated to epoch via proper motion
  Rise/set:     Standard almanac horizon corrections:
                  Sun     -0.8333° (refraction + disc)
                  Moon    +0.125°  (convention)
                  Planets -0.5667° (refraction only)
  Retrograde:   ±6h finite difference on ecliptic longitude
  Eclipses:     Topocentric sep < 0.3° → classified total/annular
  Conjunction:  6h sampling across 366-day window + 60-iteration ternary refinement
