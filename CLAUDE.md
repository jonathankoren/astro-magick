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
python3 astro_magick.py --astronomical
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
  --nighttransit  Next night-time peak for each planet + Moon rows
  --astronomical  IAU constellation positions table
  --behenian      Behenian fixed star aspects; with --nighttransit adds star peak table
  --eclipses      Next partial and total lunar/solar eclipses (up to 2053)
  --json          All requested reports as a single JSON object to stdout

All status/error messages go to stderr; stdout is clean for piping.

---

## ASTROLOGICAL POSITIONS table

Printed once per --days day. Columns:

  BODY        Planet name
  SIGN        Tropical zodiac sign
  DEG         Degrees within sign (0.0–29.9)
  Rx          " Rx" if retrograde, else blank
  RISE        Local time of horizon rise
  TRANSIT     Upper culmination. [D] = Sun up, [N] = Sun down
  SET         Local time of horizon set
  ANTITRANSIT Lower culmination (nadir). [D]/[N] as above
  NVIS        Start of night visibility window
  END         End of night visibility window

Notes:
- [D]/[N] tags are always 3 ASCII chars — no emoji
- Circumpolar bodies show `----` for RISE/SET
- All column content is pure ASCII; planet glyphs are non-padded decorative prefixes

---

## ASTRONOMICAL POSITIONS table  (--astronomical)

Second table after the astrological one, same rise/transit/set columns but:

  CONSTELLATION   IAU constellation name + % through it (e.g. `Gemini 72%`)

The ecliptic crosses 13 constellations (Ophiuchus appears between Scorpius and
Sagittarius at 247.6–266.2°). The percentage is 0 at the constellation's start
boundary and 100 at its end.

### IAU ecliptic boundary data

Source: David Asher, Human Orrery project (cantab.net), J2000 ecliptic longitudes:

  Pisces        0.000
  Aries        28.687
  Taurus       53.417
  Gemini       90.140
  Cancer      117.988
  Leo         138.038
  Virgo       173.851
  Libra       217.810
  Scorpius    241.047
  Ophiuchus   247.638
  Sagittarius 266.238
  Capricornus 299.656
  Aquarius    327.488
  Pisces      351.650  (resumes)

---

## NEXT NIGHT PEAK table  (--nighttransit)

For each body, finds the next night it is visible and shows the time of its
highest point during darkness.

### Moon rows (3 rows total)

  Moon       Next night any phase is above the horizon
  Full Moon  Next Full Moon night (elongation 165–195°)
  New Moon   Next New Moon night (elongation 345–15°, wrapping window)

All three use `find_next_night_peak_moon_phase()` with the appropriate
elongation range.

### Other bodies

  Sun               Always daytime — shown as `---`
  Mercury, Venus    30-min sampling across night window (inner planet; no night transit)
  Mars, Jupiter,    `find_next_night_transit_outer()`: enumerates all transits over
  Saturn            366 days in one Skyfield call; returns first inside a night window.
                    Handles solar conjunction gaps correctly.

### Type column

  transit     Meridian crossing falls within the night window
  night peak  Transit is in daytime; shown value is altitude max during darkness

---

## ASTRONOMICAL NEXT NIGHT PEAK table  (--astronomical --nighttransit)

Same as the planet night peak table but CONSTELLATION replaces SIGN/DEG.
Includes all three Moon rows. Printed after the astrological night peak table.

---

## BEHENIAN FIXED STARS  (--behenian)

### Daily aspects table

Shows every planet within BEHENIAN_WARN (20°) of a Behenian star.

  PLANET    Body name
  STAR      Behenian star name
  SEP       Ecliptic separation in degrees
  STATUS    "IN ORB  <6 deg" if within the traditional 6° orb, else "approaching"
  APPROACH  Skyfield-computed datetime of minimum separation (actual conjunction,
            not an estimate). "none in 1 yr" if no conjunction found within 366 days.
  VISIBLE   "night" if star is >15° from Sun (heliacal rule); "solar glare" otherwise

### Star catalog

All 15 Behenian stars are identified by Hipparcos catalog number and loaded via
`skyfield.data.hipparcos` + `skyfield.api.Star`. No hand-coded positions.
Positions are computed at runtime for the actual epoch.

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

The catalog dataframe is cached in `_HIP_DF_CACHE` / `_BEHENIAN_STAR_CACHE`
(module-level) so it is only loaded once per run.

### Conjunction time calculation

`find_conjunction_time()` uses `almanac.find_discrete` with 6-hour steps to
detect when the planet enters a <3° separation window, then refines to the
actual minimum with a 60-iteration ternary search. Retrograde planets work
correctly because Skyfield computes real motion.

---

## BEHENIAN STARS NIGHT PEAK table  (--behenian --nighttransit)

For each of the 15 Behenian stars:

  DATE               Next date the star is above the horizon at night
  TIME               Local time of its highest point during that night
  TYPE               "transit" if meridian crossing falls in night window;
                     "night peak" if transit is in daytime
  NEXT NIGHT TRANSIT Next date+time the star's meridian crossing falls within
                     a night window (scans up to 2 years forward)
  CONJ               Any planet within 6° of the star at the peak moment

---

## --eclipses table

Scans up to 2053-10-01 (DE421 limit) for:

  Lunar  partial   Next penumbral or partial lunar eclipse
  Lunar  total     Next total lunar eclipse
  Solar  partial   Next partial solar eclipse visible at location
  Solar  total     Next total or annular solar eclipse at location

Lunar eclipses: `skyfield.eclipselib.lunar_eclipses()`.
Solar eclipses: topocentric Sun–Moon separation scan + ternary search refinement.
All ephemeris calls in the solar search are wrapped in try/except for DE421 edges.

---

## --json output

All requested reports collected into one JSON object printed to stdout.
Status/error messages go to stderr so stdout is always valid JSON.

Top-level keys:

  "almanac"                    always; one entry per --days
  "astronomical"               --astronomical; one entry per --days
  "night_transit"              --nighttransit
  "astronomical_night_transit" --astronomical + --nighttransit
  "behenian"                   --behenian; one entry per --days
  "behenian_night"             --behenian + --nighttransit
  "eclipses"                   --eclipses

---

## Alignment

All column-critical content is pure ASCII. Column widths are defined in dicts
(W, NW, AW, ANW, BW, SNW) and shared between the header and every data row via
a single fmt_* closure, making misalignment structurally impossible. Planet and
zodiac glyphs appear only as non-padded decorative row prefixes.

---

## Key functions

  astronomical_constellation()       IAU boundary lookup → (name, pct_through)
  behenian_aspects()                 Planet→star aspects within BEHENIAN_WARN
  behenian_ecl_lon()                 Ecliptic lon of a Star object at time t
  build_astronomical_json()          JSON for --astronomical
  build_astronomical_night_json()    JSON for --astronomical --nighttransit
  build_behenian_json()              JSON for --behenian daily
  build_behenian_night_json()        JSON for --behenian --nighttransit
  build_daily_json()                 JSON for --almanac
  build_eclipses_json()              JSON for --eclipses
  build_nighttransit_json()          JSON for --nighttransit
  check_conjunctions()               Planet pairs within threshold degrees
  check_eclipses()                   Same-day eclipse warning for daily header
  ecl_separation()                   Shortest angular distance (0–180°)
  find_conjunction_time()            Skyfield-based actual conjunction datetime
  find_eclipses()                    Scan for lunar + solar eclipses
  find_events()                      Rise, upper transit, set, lower transit
  find_next_night_peak()             Night peak dispatcher for planets
  find_next_night_peak_moon_phase()  Moon night peak for arbitrary elongation range
  find_next_night_transit_outer()    Bulk transit scan for Mars/Jupiter/Saturn
  find_next_star_night_peak()        Night peak for a Behenian Star object
  find_next_star_night_transit()     Next nighttime transit for a Star object
  find_star_night_peak()             Highest point within one night window
  find_star_transit()                Meridian transit within a time window
  get_behenian_stars()               Load/cache Hipparcos Star objects (15 stars)
  get_body_position()                Ecliptic lon, RA, Dec at given time
  get_moon_phase()                   Phase name + illumination % from elongation
  is_retrograde()                    Detect retrograde by sampling ecl lon ±6h
  load_config()                      Read astro_magick.cfg, create if missing
  load_ephemeris()                   Load DE421, download if missing
  night_peak_altitude()              30-min sampling + ternary refinement
  night_window()                     Compute sunset → next-sunrise window
  output_json()                      Collect all reports into one JSON object
  planet_night_visibility()          Clip planet rise-set to night window
  print_astronomical_night_report()  --astronomical --nighttransit table
  print_astronomical_report()        --astronomical table
  print_behenian_night_report()      --behenian --nighttransit star peak table
  print_behenian_report()            --behenian daily aspects table
  print_eclipse_report()             --eclipses table
  print_night_transit_report()       --nighttransit planet table
  print_report()                     Daily astrological almanac table
  zodiac_sign()                      Tropical sign from ecliptic longitude
  _get_night_window_for_day()        Sunset→sunrise helper for one calendar day
  _utc_in_night()                    Check if UTC datetime falls in night window

---

## Accuracy

  Positions:    ~arcsecond level (DE421 + Skyfield apparent place)
  Fixed stars:  Hipparcos catalog positions propagated to date via proper motion
  Rise/set:     Standard almanac horizon corrections:
                  Sun     -0.8333° (refraction + disc)
                  Moon    +0.125°  (convention)
                  Planets -0.5667° (refraction only)
  Retrograde:   ±6h finite difference on ecliptic longitude
  Eclipses:     Topocentric sep < 0.3° → classified total/annular
  Conjunction:  find_discrete 6h steps + 60-iteration ternary refinement
