# CLAUDE.md — Astro Magick

Developer notes for `astro_magick.py`.

---

## Overview

Command-line astronomical almanac for magickal practice.
Powered by **Skyfield + JPL DE421 ephemeris** (~17 MB, downloaded on first run).

```
python3 astro_magick.py
python3 astro_magick.py --date 2026-06-21
python3 astro_magick.py --lat 51.5074 --lon -0.1278 --timezone Europe/London
python3 astro_magick.py --days 7
python3 astro_magick.py --nighttransit
python3 astro_magick.py --eclipses
```

---

## Dependencies

```
pip install skyfield
```

DE421 covers 1900-2050. For dates outside that range, switch to `de440.bsp`
in `load_ephemeris()`.

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

## Command-line options

  --lat          Latitude (decimal degrees, + = North)
  --lon          Longitude (decimal degrees, + = East)
  --timezone     IANA timezone name
  --date         Start date YYYY-MM-DD (default: today)
  --days N       Number of consecutive daily reports
  --nighttransit Show next night-time peak altitude for each planet
  --eclipses     Show next partial and total lunar/solar eclipses

---

## Main daily table columns

  sym  BODY       SIGN           DEG    Rx  RISE  TRANSIT   SET  ANTITRANSIT  NVIS   END

  sym         Decorative planet glyph (outside column, not padded)
  BODY        Planet name
  SIGN        Zodiac sign name (pure ASCII -- no glyph, for alignment)
  DEG         Degrees within sign (0.0-29.9)
  Rx          " Rx" if planet is retrograde, else blank
  RISE        Time of horizon rise
  TRANSIT     Upper culmination (highest point). [D]=Sun up, [N]=Sun down
  SET         Time of horizon set
  ANTITRANSIT Lower culmination (nadir). [D]/[N] as above
  NVIS        Start of night visibility window (rise-set overlap with darkness)
  END         End of night visibility window

Notes:
- [D] / [N] tags are always 3 ASCII chars -- no emoji, no ambiguous-width chars
- Circumpolar bodies show ---- for RISE/SET; night vis still computed
- All column content is pure ASCII; decorative glyphs are fixed non-padded prefixes

---

## --nighttransit table

For each planet, finds the next night when it reaches its peak altitude:

  Moon              Scans until the next Full Moon (elongation 165-195 deg)
                    that occurs at night
  Mars/Jupiter/Saturn  find_next_night_transit_outer(): enumerates all transits
                    over 366 days in one Skyfield call, returns first one in a
                    night window. Handles conjunction gaps correctly.
  Mercury/Venus     Inner planets never transit at night; finds next night-time
                    altitude peak (dusk/dawn apparition) via 30-min sampling
  Sun               Always daytime -- shown as ---

Type column:
  transit     Body crosses meridian at highest point during darkness
  night peak  Transit is in daytime; this is altitude max during darkness

---

## --eclipses table

Scans up to 2 years forward for:

  Lunar  partial   Next penumbral or partial lunar eclipse
  Lunar  total     Next total lunar eclipse
  Solar  partial   Next partial solar eclipse visible at location
  Solar  total     Next total or annular solar eclipse at location

Lunar eclipse times = geocentric peak (globally visible if Moon is up).
Solar eclipse times = moment of minimum topocentric Sun-Moon separation.

---

## Alignment

All column-critical content is pure ASCII because:
- Zodiac glyphs (Aries-Pisces) are Unicode Wide (W) = always 2 terminal columns
- Planet glyphs are Ambiguous (A) = 1 or 2 columns depending on terminal
- Emoji have inconsistent multi-column rendering

[D]/[N] replaced emoji for transit indicators.
" Rx" (3 ASCII chars) replaced the Rx glyph for retrograde.
Planet/zodiac symbols appear only as decorative non-padded row prefixes.

---

## Key functions

  load_config()                    Read astro_magick.cfg, create if missing
  load_ephemeris()                 Load DE421, download if missing
  get_body_position()              Ecliptic lon, RA, Dec at given time
  is_retrograde()                  Detect retrograde by sampling ecl lon +/-6h
  find_events()                    Rise, upper transit, set, lower transit
  get_moon_phase()                 Phase name + illumination % from elongation
  check_eclipses()                 Same-day eclipse warnings for daily header
  check_conjunctions()             Planet pairs within threshold degrees
  night_window()                   Compute sunset -> next-sunrise window
  planet_night_visibility()        Clip planet rise-set to night window
  night_peak_altitude()            30-min sampling + ternary refinement
  find_next_night_transit_outer()  Bulk transit scan for Mars/Jupiter/Saturn
  find_next_night_peak()           Dispatcher to above or sampling
  find_eclipses()                  Scan 2 years for lunar + solar eclipses
  print_report()                   Daily almanac table
  print_night_transit_report()     Next night peak table
  print_eclipse_report()           Eclipse forecast table

---

## Accuracy

  Positions:  ~arcsecond level (DE421 + Skyfield apparent place)
  Rise/set:   Standard almanac horizon corrections:
                Sun    -0.8333 deg (refraction + disc)
                Moon   +0.125  deg (convention)
                Planets -0.5667 deg (refraction only)
  Retrograde: +/-6h finite difference on ecliptic longitude
  Eclipse:    Topocentric Sep < 0.3 deg -> classified total/annular

---

## Zodiac

Tropical only (ecliptic longitude from J2000 frame). No sidereal support.
