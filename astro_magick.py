#!/usr/bin/env python3
"""
✦ ASTRO MAGICK ✦
Command-line astronomical almanac for magickal practice.
Powered by Skyfield + JPL ephemeris DE421 for high-accuracy positions.

On first run, skyfield will download de421.bsp (~17MB) to the current directory.

Usage:
    python3 astro_magick.py
    python3 astro_magick.py --date 2026-03-21
    python3 astro_magick.py --lat 51.5074 --lon -0.1278 --timezone Europe/London
    python3 astro_magick.py --days 7
"""

import argparse
import configparser
import datetime
import json
import math
import pathlib
import sys
from zoneinfo import ZoneInfo

try:
    from skyfield import almanac
    from skyfield.api import load, wgs84, N, E
    from skyfield.framelib import ecliptic_frame
    from skyfield.eclipselib import lunar_eclipses
except ImportError:
    print("ERROR: skyfield is not installed.", file=sys.stderr)
    print("Install it with:  pip install skyfield")
    sys.exit(1)

# ─────────────────────────────────────────────
#  CONSTANTS
# ─────────────────────────────────────────────

ZODIAC = [
    ("Aries",        0,  30), ("Taurus",      30,  60), ("Gemini",      60,  90),
    ("Cancer",      90, 120), ("Leo",        120, 150), ("Virgo",      150, 180),
    ("Libra",      180, 210), ("Scorpio",    210, 240), ("Sagittarius",240, 270),
    ("Capricorn",  270, 300), ("Aquarius",   300, 330), ("Pisces",     330, 360),
]

ZODIAC_SYMBOLS = {
    "Aries": "♈", "Taurus": "♉", "Gemini": "♊", "Cancer": "♋",
    "Leo": "♌", "Virgo": "♍", "Libra": "♎", "Scorpio": "♏",
    "Sagittarius": "♐", "Capricorn": "♑", "Aquarius": "♒", "Pisces": "♓",
}

PLANET_SYMBOLS = {
    "Sun": "☉", "Moon": "☽", "Mercury": "☿", "Venus": "♀",
    "Mars": "♂", "Jupiter": "♃", "Saturn": "♄",
}

SKYFIELD_NAMES = {
    "Sun":     "sun",
    "Moon":    "moon",
    "Mercury": "mercury",
    "Venus":   "venus",
    "Mars":    "mars",
    "Jupiter": "jupiter barycenter",
    "Saturn":  "saturn barycenter",
}

# ─────────────────────────────────────────────
#  HELPERS
# ─────────────────────────────────────────────

def zodiac_sign(lon_deg):
    lon = lon_deg % 360.0
    for name, start, end in ZODIAC:
        if start <= lon < end:
            return name
    raise ValueError(
        f"Could not determine zodiac sign for ecliptic longitude "
        f"{lon_deg:.6f}° (normalized: {lon:.6f}°)"
    )

def fmt_time(dt_utc, tz):
    if dt_utc is None:
        return "---"
    return dt_utc.astimezone(tz).strftime("%H:%M")


# ─────────────────────────────────────────────
#  EPHEMERIS SETUP
# ─────────────────────────────────────────────

def load_ephemeris():
    try:
        eph = load("de421.bsp")
    except Exception as e:
        print(f"ERROR loading ephemeris: {e}", file=sys.stderr)
        print("Make sure you have an internet connection for the first run.")
        sys.exit(1)
    ts = load.timescale()
    return eph, ts

# ─────────────────────────────────────────────
#  BODY POSITION
# ─────────────────────────────────────────────

def get_body_position(eph, ts, name, t, observer):
    """Returns (ecl_lon_deg, ra_deg, dec_deg) apparent from observer at time t."""
    body = eph[SKYFIELD_NAMES[name]]
    astrometric = observer.at(t).observe(body).apparent()
    _, ecl_lon, _ = astrometric.frame_latlon(ecliptic_frame)
    ra, dec, _ = astrometric.radec()
    return ecl_lon.degrees % 360.0, ra._degrees % 360.0, dec.degrees

# ─────────────────────────────────────────────
#  RETROGRADE DETECTION
# ─────────────────────────────────────────────

# The Sun and Moon never go retrograde; no need to check them.
RETROGRADE_BODIES = {"Mercury", "Venus", "Mars", "Jupiter", "Saturn"}

def is_retrograde(eph, ts, name, t):
    """
    Return True if the body is currently in retrograde motion.

    Retrograde means the planet's apparent geocentric ecliptic longitude is
    decreasing (moving westward against the stars) rather than increasing.
    We detect this by comparing longitude 6 hours before and after t and
    checking the sign of the difference, accounting for the 0°/360° wrap.
    """
    if name not in RETROGRADE_BODIES:
        return False

    body   = eph[SKYFIELD_NAMES[name]]
    earth  = eph["earth"]
    dt     = 6 / 24  # 6 hours in days

    def ecl_lon_at(jd_offset):
        t_sample = ts.tt_jd(t.tt + jd_offset)
        astr = earth.at(t_sample).observe(body).apparent()
        _, lon, _ = astr.frame_latlon(ecliptic_frame)
        return lon.degrees

    lon_before = ecl_lon_at(-dt)
    lon_after  = ecl_lon_at(+dt)

    # Unwrap the difference across the 0°/360° boundary
    diff = (lon_after - lon_before + 180) % 360 - 180
    return diff < 0


# ─────────────────────────────────────────────
#  RISE / CULMINATION / SET  (via skyfield almanac)
# ─────────────────────────────────────────────

def find_events(eph, ts, name, date, observer, topos, tz):
    """
    Find rise, upper culmination (celestial noon), set, and lower culmination
    (celestial midnight) for a body on a given local date.

    observer = eph["earth"] + topos  (the full geocentric + surface vector)
    topos    = wgs84.latlon(...)     (the bare surface position)

    Upper culmination  — body crosses the meridian at its highest point above
                         the horizon. This is the body's personal "noon".
    Lower culmination  — body crosses the meridian at its lowest point below
                         the horizon (nadir). This is the body's personal "midnight".

    Returns (rise_utc, culmination_utc, set_utc, nadir_utc, status).
    Status is 'normal', 'circumpolar', or 'never'.
    """
    local_start = datetime.datetime(date.year, date.month, date.day,  0, 0, 0, tzinfo=tz)
    local_end   = datetime.datetime(date.year, date.month, date.day, 23,59,59, tzinfo=tz)
    t0 = ts.from_datetime(local_start.astimezone(datetime.timezone.utc))
    t1 = ts.from_datetime(local_end.astimezone(datetime.timezone.utc))

    body = eph[SKYFIELD_NAMES[name]]

    # Horizon thresholds matching standard almanac conventions:
    #   Sun:     -0.8333° (refraction + solar disc upper limb)
    #   Moon:     0.125° (refraction + lunar disc, standard convention)
    #   Planets: -0.5667° (refraction only, point source)
    if name == "Sun":
        horizon_deg = -0.8333
    elif name == "Moon":
        horizon_deg = 0.125
    else:
        horizon_deg = -0.5667

    # --- Rise and Set ---
    # find_risings/find_settings take the combined observer vector and return
    # (times, flags) where flag=False means the body never actually crossed the
    # horizon (circumpolar or never-rises).
    rise_times, rise_flags = almanac.find_risings(observer, body, t0, t1,
                                                   horizon_degrees=horizon_deg)
    set_times,  set_flags  = almanac.find_settings(observer, body, t0, t1,
                                                    horizon_degrees=horizon_deg)

    # Discard grazing events where the body didn't truly cross the horizon.
    real_rises = [t for t, f in zip(rise_times, rise_flags) if f]
    real_sets  = [t for t, f in zip(set_times,  set_flags)  if f]

    rise_utc = real_rises[0].utc_datetime() if real_rises else None
    set_utc  = real_sets[0].utc_datetime()  if real_sets  else None

    # --- Upper culmination (celestial noon) and lower culmination (celestial midnight) ---
    # meridian_transits returns event=1 for upper transit, event=0 for lower transit.
    f_meridian = almanac.meridian_transits(eph, body, topos)
    trans_times, trans_events = almanac.find_discrete(t0, t1, f_meridian)

    culmination_utc = None  # highest point — celestial noon
    nadir_utc       = None  # lowest point  — celestial midnight
    for t, event in zip(trans_times, trans_events):
        if event == 1 and culmination_utc is None:
            culmination_utc = t.utc_datetime()
        elif event == 0 and nadir_utc is None:
            nadir_utc = t.utc_datetime()

    # --- Circumpolar / never-rises detection ---
    if rise_utc is None and set_utc is None:
        local_noon = datetime.datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tz)
        t_noon = ts.from_datetime(local_noon.astimezone(datetime.timezone.utc))
        alt, _, _ = observer.at(t_noon).observe(body).apparent().altaz()
        status = "circumpolar" if alt.degrees > horizon_deg else "never"
    else:
        status = "normal"

    return rise_utc, culmination_utc, set_utc, nadir_utc, status

# ─────────────────────────────────────────────
#  MOON PHASE
# ─────────────────────────────────────────────

def get_moon_phase(eph, ts, t):
    earth = eph["earth"]
    _, sun_lon, _  = earth.at(t).observe(eph["sun"]).apparent().frame_latlon(ecliptic_frame)
    _, moon_lon, _ = earth.at(t).observe(eph["moon"]).apparent().frame_latlon(ecliptic_frame)
    elongation = (moon_lon.degrees - sun_lon.degrees) % 360.0
    illumination = (1 - math.cos(math.radians(elongation))) / 2 * 100

    if elongation < 22.5:    phase = "🌑 New Moon"
    elif elongation < 67.5:  phase = "🌒 Waxing Crescent"
    elif elongation < 112.5: phase = "🌓 First Quarter"
    elif elongation < 157.5: phase = "🌔 Waxing Gibbous"
    elif elongation < 202.5: phase = "🌕 Full Moon"
    elif elongation < 247.5: phase = "🌖 Waning Gibbous"
    elif elongation < 292.5: phase = "🌗 Last Quarter"
    elif elongation < 337.5: phase = "🌘 Waning Crescent"
    else:                    phase = "🌑 New Moon"

    return phase, round(illumination, 1)

# ─────────────────────────────────────────────
#  ECLIPSE DETECTION  (via skyfield almanac)
# ─────────────────────────────────────────────

LUNAR_ECLIPSE_TYPES = {0: "Penumbral", 1: "Partial", 2: "Total"}

def check_eclipses(eph, ts, date):
    """
    Check for any solar or lunar eclipses on the given date.
    Lunar eclipses use skyfield.eclipselib.lunar_eclipses().
    Solar eclipses use almanac.find_discrete() with a custom shadow function,
    since skyfield has no built-in solar eclipse finder.
    Returns a list of warning strings (empty if none).
    """
    warnings = []

    # Search ±2 days so eclipses peaking near midnight aren't missed.
    d0 = date - datetime.timedelta(days=2)
    d1 = date + datetime.timedelta(days=2)
    t0 = ts.utc(d0.year, d0.month, d0.day)
    t1 = ts.utc(d1.year, d1.month, d1.day)

    # --- Lunar eclipses (skyfield.eclipselib) ---
    eclipse_times, eclipse_types, _ = lunar_eclipses(t0, t1, eph)
    for t_ecl, etype in zip(eclipse_times, eclipse_types):
        if t_ecl.utc_datetime().date() == date:
            kind = LUNAR_ECLIPSE_TYPES.get(int(etype), "Unknown")
            warnings.append(f"🌕 LUNAR ECLIPSE today — {kind}")

    # --- Solar eclipses (almanac.find_discrete with custom function) ---
    # A solar eclipse occurs when the Moon's penumbral or umbral shadow crosses
    # Earth's surface, i.e. when the Sun–Moon angular separation (as seen from
    # Earth's centre) is less than the sum of their apparent radii (~1.5°).
    # We build a discrete step function that is 1 inside that threshold and 0
    # outside, then let find_discrete locate the entry/exit crossings.
    earth = eph["earth"]
    sun   = eph["sun"]
    moon  = eph["moon"]

    def solar_eclipse_possible(t):
        e = earth.at(t)
        s = e.observe(sun).apparent()
        m = e.observe(moon).apparent()
        return m.separation_from(s).degrees < 1.5

    solar_eclipse_possible.step_days = 1.0  # check at most once per day

    eclipse_events, eclipse_flags = almanac.find_discrete(t0, t1, solar_eclipse_possible)
    for t_ecl, flag in zip(eclipse_events, eclipse_flags):
        # flag=1 means we just entered the eclipse window (separation dropped below 1.5°)
        if flag and t_ecl.utc_datetime().date() == date:
            warnings.append("☀️  SOLAR ECLIPSE possible today — Moon's shadow near Earth")

    return warnings

# ─────────────────────────────────────────────
#  CONJUNCTIONS
# ─────────────────────────────────────────────

def check_conjunctions(positions, threshold=8.0):
    names = list(positions.keys())
    results = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            diff = abs(positions[a] - positions[b])
            if diff > 180:
                diff = 360 - diff
            if diff <= threshold:
                results.append((a, b, round(diff, 2)))
    return results

# ─────────────────────────────────────────────
#  NIGHT VISIBILITY & NEXT-NIGHT-TRANSIT HELPERS
# ─────────────────────────────────────────────

def night_window(sunrise_utc, sunset_utc, date, tz, ts, observer, sun_body):
    """
    Return (night_start_utc, night_end_utc) for the night that BEGINS on `date`
    (i.e. after that day's sunset until the following morning's sunrise).

    sunset_utc  — today's sunset (may be None for polar regions)
    sunrise_utc — today's sunrise (may be None); we need tomorrow's sunrise

    Returns (None, None) for polar summer (no night) or polar winter (no day).
    """
    # Compute tomorrow's sunrise
    tomorrow = date + datetime.timedelta(days=1)
    t_tmrw_start = ts.from_datetime(
        datetime.datetime(tomorrow.year, tomorrow.month, tomorrow.day,
                          0, 0, 0, tzinfo=tz).astimezone(datetime.timezone.utc))
    t_tmrw_end = ts.from_datetime(
        datetime.datetime(tomorrow.year, tomorrow.month, tomorrow.day,
                          23, 59, 59, tzinfo=tz).astimezone(datetime.timezone.utc))
    tmrw_rise_times, tmrw_rise_flags = almanac.find_risings(
        observer, sun_body, t_tmrw_start, t_tmrw_end, horizon_degrees=-0.8333)
    real_tmrw_rises = [t.utc_datetime() for t, f in zip(tmrw_rise_times, tmrw_rise_flags) if f]
    tomorrow_sunrise_utc = real_tmrw_rises[0] if real_tmrw_rises else None

    night_start = sunset_utc          # tonight's sunset
    night_end   = tomorrow_sunrise_utc # tomorrow's sunrise

    if night_start is None or night_end is None:
        return None, None
    if night_end <= night_start:
        return None, None
    return night_start, night_end


def planet_night_visibility(rise_utc, set_utc, status, night_start, night_end):
    """
    Given a planet's rise/set and the night window, return
    (vis_start_utc, vis_end_utc) for the overlap, or (None, None) if none.
    """
    if night_start is None or night_end is None:
        return None, None
    if status == "never":
        return None, None
    if status == "circumpolar":
        # Always up — visible the whole night
        return night_start, night_end

    # Normal case: clip [rise, set] to [night_start, night_end]
    if rise_utc is None or set_utc is None:
        return None, None
    vis_start = max(rise_utc, night_start)
    vis_end   = min(set_utc,  night_end)
    if vis_end <= vis_start:
        return None, None
    return vis_start, vis_end


def night_peak_altitude(ts, observer, body, night_start_utc, night_end_utc):
    """
    Find the time of peak altitude within the night window [night_start_utc, night_end_utc]
    using a ternary search. Returns the skyfield Time object of peak altitude,
    or None if the night window is invalid.
    """
    if night_start_utc is None or night_end_utc is None:
        return None
    t_lo = ts.from_datetime(night_start_utc.replace(tzinfo=datetime.timezone.utc))
    t_hi = ts.from_datetime(night_end_utc.replace(tzinfo=datetime.timezone.utc))
    for _ in range(60):
        span = t_hi.tt - t_lo.tt
        if span < 1e-8:
            break
        m1 = ts.tt_jd(t_lo.tt + span / 3)
        m2 = ts.tt_jd(t_hi.tt - span / 3)
        a1 = observer.at(m1).observe(body).apparent().altaz()[0].degrees
        a2 = observer.at(m2).observe(body).apparent().altaz()[0].degrees
        if a1 < a2:
            t_lo = m1
        else:
            t_hi = m2
    return ts.tt_jd((t_lo.tt + t_hi.tt) / 2)


def find_next_night_transit_outer(eph, ts, name, start_date, observer, topos, tz, max_days=366):
    """
    For slow outer planets (Jupiter, Saturn) whose transit drifts only ~4
    minutes per day: collect ALL upper transits over max_days in a single
    skyfield call, then walk through them in order and return the first one
    whose UTC time falls within the local night window (sunset → next sunrise).

    This avoids the day-by-day loop that breaks when the planet is near solar
    conjunction and below the horizon every night for weeks at a time.
    """
    earth    = eph["earth"]
    sun_body = eph[SKYFIELD_NAMES["Sun"]]
    body     = eph[SKYFIELD_NAMES[name]]

    t_start = ts.from_datetime(
        datetime.datetime(start_date.year, start_date.month, start_date.day,
                          0, 0, 0, tzinfo=tz).astimezone(datetime.timezone.utc))
    t_end   = ts.from_datetime(
        datetime.datetime(start_date.year, start_date.month, start_date.day,
                          0, 0, 0, tzinfo=tz).astimezone(datetime.timezone.utc)
        + datetime.timedelta(days=max_days))

    all_transits = almanac.find_transits(observer, body, t_start, t_end)

    for transit_t in all_transits:
        transit_utc = transit_t.utc_datetime()
        # The local date on which this transit falls
        local_dt = transit_utc.astimezone(tz)
        # Night window starts the evening of the day before the transit
        # if the transit is after midnight, or the same evening if before.
        # Simplest: check the night whose evening is transit_date - 1 day
        # AND the night whose evening is transit_date itself.
        for day_offset in (-1, 0):
            d = local_dt.date() + datetime.timedelta(days=day_offset)
            local_start = datetime.datetime(d.year, d.month, d.day, 0, 0, 0, tzinfo=tz)
            local_end   = datetime.datetime(d.year, d.month, d.day, 23, 59, 59, tzinfo=tz)
            t0 = ts.from_datetime(local_start.astimezone(datetime.timezone.utc))
            t1 = ts.from_datetime(local_end.astimezone(datetime.timezone.utc))

            ss_times, ss_flags = almanac.find_settings(observer, sun_body, t0, t1,
                                                        horizon_degrees=-0.8333)
            sr_times, sr_flags = almanac.find_risings(observer, sun_body, t0, t1,
                                                       horizon_degrees=-0.8333)
            real_ss = [t.utc_datetime() for t, f in zip(ss_times, ss_flags) if f]
            real_sr = [t.utc_datetime() for t, f in zip(sr_times, sr_flags) if f]
            sunset_utc  = real_ss[0] if real_ss else None
            sunrise_utc = real_sr[0] if real_sr else None

            night_start_utc, night_end_utc = night_window(
                sunrise_utc, sunset_utc, d, tz, ts, observer, sun_body)

            if night_start_utc is None or night_end_utc is None:
                continue
            if night_start_utc <= transit_utc <= night_end_utc:
                # Found a night transit — collect position and return
                astr = earth.at(transit_t).observe(body).apparent()
                _, ecl_lon_obj, _ = astr.frame_latlon(ecliptic_frame)
                ecl_lon = ecl_lon_obj.degrees % 360.0
                return dict(date=d, peak_utc=transit_utc, ecl_lon=ecl_lon,
                            is_transit=True, full_moon=False)

    return None


def find_next_night_peak(eph, ts, name, start_date, observer, topos, tz, max_days=366):
    """
    Scan forward from start_date to find the next night the body is visible,
    and return its peak altitude time during that night window.

    The night window is sunset(day d) → sunrise(day d+1).  We search for an
    upper transit *within that exact window* (not within the calendar day) so
    that post-midnight transits are found correctly.

    If a transit falls in the night window  → is_transit=True, peak=transit.
    If no transit in the window             → sample altitude at 30-minute
        intervals across the night window and take the highest above-horizon
        sample.  This handles planets that are not unimodal in altitude (e.g.
        a planet that sets just after sunset and rises again before dawn —
        ternary search would converge on the wrong hump).

    For the Moon the function scans until it finds a Full Moon night
    (elongation 165–195°).  For all other bodies it returns the first night
    the body has any above-horizon time.

    Returns a dict:
        date        — local calendar date on which the night begins
        peak_utc    — UTC datetime of peak altitude during the night
        ecl_lon     — ecliptic longitude at peak (degrees)
        is_transit  — True if peak is a meridian transit
        full_moon   — True if Moon is near Full at peak time
    Returns None if no qualifying night is found within max_days.
    """
    # Outer planets (Mars, Jupiter, Saturn) can oppose the Sun and transit at
    # midnight.  Near conjunction they're invisible for weeks.  Enumerate all
    # their transits over the full window in one skyfield call rather than the
    # day-by-day loop, which breaks near conjunction.
    if name in ("Mars", "Jupiter", "Saturn"):
        return find_next_night_transit_outer(
            eph, ts, name, start_date, observer, topos, tz, max_days)

    earth    = eph["earth"]
    sun_body = eph[SKYFIELD_NAMES["Sun"]]
    body     = eph[SKYFIELD_NAMES[name]]

    for offset in range(max_days):
        d = start_date + datetime.timedelta(days=offset)

        # Build the day window to find sunset
        local_start = datetime.datetime(d.year, d.month, d.day, 0, 0, 0, tzinfo=tz)
        local_end   = datetime.datetime(d.year, d.month, d.day, 23, 59, 59, tzinfo=tz)
        t0 = ts.from_datetime(local_start.astimezone(datetime.timezone.utc))
        t1 = ts.from_datetime(local_end.astimezone(datetime.timezone.utc))

        ss_times, ss_flags = almanac.find_settings(observer, sun_body, t0, t1,
                                                    horizon_degrees=-0.8333)
        real_ss = [t.utc_datetime() for t, f in zip(ss_times, ss_flags) if f]
        sunset_utc = real_ss[0] if real_ss else None

        sr_times, sr_flags = almanac.find_risings(observer, sun_body, t0, t1,
                                                   horizon_degrees=-0.8333)
        real_sr = [t.utc_datetime() for t, f in zip(sr_times, sr_flags) if f]
        sunrise_utc = real_sr[0] if real_sr else None

        night_start_utc, night_end_utc = night_window(
            sunrise_utc, sunset_utc, d, tz, ts, observer, sun_body)

        if night_start_utc is None:
            continue  # polar summer — no night this day

        # Search for a transit within the night window (sunset → next sunrise).
        # Using the night window as the search range catches post-midnight transits.
        t_night_start = ts.from_datetime(
            night_start_utc.replace(tzinfo=datetime.timezone.utc))
        t_night_end   = ts.from_datetime(
            night_end_utc.replace(tzinfo=datetime.timezone.utc))

        trans_times = almanac.find_transits(observer, body, t_night_start, t_night_end)
        if len(trans_times) > 0:
            peak_t     = trans_times[0]
            is_transit = True
        else:
            # No transit in the night window.  Sample altitude at 30-min
            # intervals across the night and find the highest above-horizon
            # point.  Sampling avoids the unimodal assumption of ternary search
            # — a planet near solar conjunction may set right after sunset and
            # rise again just before dawn, giving two humps in the altitude
            # curve; ternary search would converge on the trough between them.
            night_duration_days = (
                (night_end_utc - night_start_utc).total_seconds() / 86400.0)
            n_steps   = max(2, int(night_duration_days * 48))  # ~30-min steps
            step_size = night_duration_days / n_steps

            best_alt = -90.0
            best_t   = None
            for step in range(n_steps + 1):
                t_sample = ts.tt_jd(t_night_start.tt + step * step_size)
                alt_deg  = observer.at(t_sample).observe(body).apparent().altaz()[0].degrees
                if alt_deg > best_alt:
                    best_alt = alt_deg
                    best_t   = t_sample

            if best_t is None or best_alt <= 0:
                continue  # planet below horizon all night — try next day

            # Refine the best sample with a short ternary search in a ±1-hour
            # window around it, now safe because we've isolated the correct hump
            lo_tt = best_t.tt - 1 / 24
            hi_tt = best_t.tt + 1 / 24
            # Clamp to night window
            lo_tt = max(lo_tt, t_night_start.tt)
            hi_tt = min(hi_tt, t_night_end.tt)
            for _ in range(40):
                span = hi_tt - lo_tt
                if span < 1e-7:
                    break
                m1 = ts.tt_jd(lo_tt + span / 3)
                m2 = ts.tt_jd(hi_tt - span / 3)
                a1 = observer.at(m1).observe(body).apparent().altaz()[0].degrees
                a2 = observer.at(m2).observe(body).apparent().altaz()[0].degrees
                if a1 < a2:
                    lo_tt = lo_tt + span / 3
                else:
                    hi_tt = hi_tt - span / 3
            peak_t     = ts.tt_jd((lo_tt + hi_tt) / 2)
            is_transit = False

        peak_utc = peak_t.utc_datetime()

        # Ecliptic longitude at peak
        astr = earth.at(peak_t).observe(body).apparent()
        _, ecl_lon_obj, _ = astr.frame_latlon(ecliptic_frame)
        ecl_lon = ecl_lon_obj.degrees % 360.0

        # Full moon check
        full_moon = False
        if name == "Moon":
            _, sun_lon_obj,  _ = (earth.at(peak_t)
                                  .observe(sun_body).apparent()
                                  .frame_latlon(ecliptic_frame))
            elongation = (ecl_lon - sun_lon_obj.degrees) % 360.0
            full_moon  = 165 <= elongation <= 195
            if not full_moon:
                continue  # keep scanning until we find a Full Moon night

        return dict(date=d, peak_utc=peak_utc, ecl_lon=ecl_lon,
                    is_transit=is_transit, full_moon=full_moon)

    return None


# ─────────────────────────────────────────────
#  MAIN REPORT
# ─────────────────────────────────────────────

def print_report(eph, ts, date, lat, lon, tz):
    lon_sign = E if lon >= 0 else -E
    topos    = wgs84.latlon(lat * N, abs(lon) * lon_sign)
    observer = eph["earth"] + topos

    local_noon = datetime.datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tz)
    utc_offset = local_noon.utcoffset().total_seconds() / 3600.0
    t_noon     = ts.from_datetime(local_noon.astimezone(datetime.timezone.utc))

    hdr_width = 70

    print()
    print()
    print("=" * hdr_width)
    print(f"  ASTRO MAGICK ALMANAC  --  {date.strftime('%A, %B %d, %Y')}")
    print("=" * hdr_width)
    print(f"  Location: {lat:.4f}N, {lon:.4f}E  |  {tz.key}  (UTC{utc_offset:+.1f})")
    print(f"  All times shown in local time")
    print()
    print()

    phase, illum = get_moon_phase(eph, ts, t_noon)
    print(f"  MOON PHASE:  {phase}  ({illum}% illuminated)")
    print()
    print()

    eclipse_warnings = check_eclipses(eph, ts, date)
    for warning in eclipse_warnings:
        print(f"  ⚠️  {warning}")
    if eclipse_warnings:
        print()

    local_start = datetime.datetime(date.year, date.month, date.day,  0, 0, 0, tzinfo=tz)
    local_end   = datetime.datetime(date.year, date.month, date.day, 23,59,59, tzinfo=tz)
    t0 = ts.from_datetime(local_start.astimezone(datetime.timezone.utc))
    t1 = ts.from_datetime(local_end.astimezone(datetime.timezone.utc))

    # Pre-compute sun rise/set and night window for this date
    sun_body = eph[SKYFIELD_NAMES["Sun"]]
    sun_rise_times, sun_rise_flags = almanac.find_risings(observer,  sun_body, t0, t1, horizon_degrees=-0.8333)
    sun_set_times,  sun_set_flags  = almanac.find_settings(observer, sun_body, t0, t1, horizon_degrees=-0.8333)
    real_sunrises = [t.utc_datetime() for t, f in zip(sun_rise_times, sun_rise_flags) if f]
    real_sunsets  = [t.utc_datetime() for t, f in zip(sun_set_times,  sun_set_flags)  if f]
    sunrise_utc = real_sunrises[0] if real_sunrises else None
    sunset_utc  = real_sunsets[0]  if real_sunsets  else None

    night_start, night_end = night_window(
        sunrise_utc, sunset_utc, date, tz, ts, observer, sun_body)

    def sun_up_at(dt_utc):
        """True if the Sun is above the horizon at the given UTC datetime."""
        if dt_utc is None:
            return None
        if sunrise_utc is not None and sunset_utc is not None:
            return sunrise_utc <= dt_utc <= sunset_utc
        t_check = ts.from_datetime(dt_utc.replace(tzinfo=datetime.timezone.utc))
        alt, _, _ = observer.at(t_check).observe(sun_body).apparent().altaz()
        return alt.degrees > -0.8333

    def transit_tag(dt_utc):
        """[D] if Sun is up (daytime), [N] if not (nighttime). Pure ASCII, always 3 chars."""
        if dt_utc is None:
            return "   "
        return "[D]" if sun_up_at(dt_utc) else "[N]"

    # All column content is pure ASCII for reliable alignment.
    # Decorative planet symbol printed as fixed prefix, not inside padded field.
    # [D] = daytime (Sun above horizon), [N] = nighttime
    # Column widths — header and data rows use identical format spec
    # so they cannot diverge.  All content is plain ASCII.
    W = dict(body=11, sign=13, deg=6, rx=3, rise=5,
             transit=9, set_=6, anti=11, nvis=5, end=5)

    def fmt_row(body, sign, deg, rx, rise, transit, set_, anti, nvis, end, note=''):
        return (f"  {body:<{W['body']}} {sign:<{W['sign']}} {deg:>{W['deg']}}"
                f"{rx:<{W['rx']}} {rise:>{W['rise']}} {transit:>{W['transit']}}"
                f" {set_:>{W['set_']}} {anti:>{W['anti']}}  {nvis:>{W['nvis']}}  {end:>{W['end']}}"
                + note)

    hdr = fmt_row('BODY','SIGN','DEG','','RISE','TRANSIT','SET','ANTITRANSIT','NVIS','END')
    width = len(hdr)
    bar   = "-" * width
    print(bar)
    print(hdr)
    print(bar)

    bodies_lon = {}
    bodies_rx  = {}

    for name in ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]:
        ecl_lon, _, _ = get_body_position(eph, ts, name, t_noon, observer)
        sign        = zodiac_sign(ecl_lon)
        symbol      = PLANET_SYMBOLS[name]
        retrograde  = is_retrograde(eph, ts, name, t_noon)
        rx_str      = " Rx" if retrograde else "   "
        deg_in_sign = ecl_lon % 30

        rise_utc, culmination_utc, set_utc, nadir_utc, status = find_events(
            eph, ts, name, date, observer, topos, tz
        )

        if status == "circumpolar":
            rise_str = set_str = " ----"
            note = "  (*) circumpolar"
        elif status == "never":
            rise_str = set_str = " ----"
            note = "  (-) never rises"
        else:
            rise_str = fmt_time(rise_utc, tz)
            set_str  = fmt_time(set_utc,  tz)
            note = ""

        cul_tag = transit_tag(culmination_utc)
        nad_tag = transit_tag(nadir_utc)
        cul_str = f"{cul_tag}{fmt_time(culmination_utc, tz)}"   # e.g. "[D]12:21"
        nad_str = f"{nad_tag}{fmt_time(nadir_utc, tz)}"         # e.g. "[N]00:21"

        vis_start, vis_end = planet_night_visibility(
            rise_utc, set_utc, status, night_start, night_end)
        vis_start_str = fmt_time(vis_start, tz) if vis_start else " ----"
        vis_end_str   = fmt_time(vis_end,   tz) if vis_end   else " ----"

        print(fmt_row(name, sign, f"{deg_in_sign:.1f}", rx_str,
                       rise_str, cul_str, set_str, nad_str,
                       vis_start_str, vis_end_str, note))
        bodies_lon[name] = ecl_lon
        bodies_rx[name]  = retrograde

    print(bar)

    conj = check_conjunctions(bodies_lon)
    if conj:
        print()
        print("  CONJUNCTIONS & CLOSE ASPECTS (within 8 deg):")
        for a, b, diff in conj:
            line = f"    {PLANET_SYMBOLS[a]} {a}  conj  {PLANET_SYMBOLS[b]} {b}  -- {diff} deg apart"
            if diff <= 3.0:
                line += "  [very tight]"
            print(line)

    print()
    print()
    print("=" * hdr_width)
    print()
    print()

# ─────────────────────────────────────────────
#  NEXT-NIGHT-TRANSIT REPORT
# ─────────────────────────────────────────────

def print_night_transit_report(eph, ts, start_date, lat, lon, tz):
    """
    For each planet, find the next night on which it is visible and show the
    time of its highest point in the sky during that night.

    If the planet's upper transit falls within the night window, that is the
    peak — labelled TRANSIT. Otherwise the peak is wherever the planet reaches
    maximum altitude during darkness — labelled NIGHT PEAK (not a transit).
    """
    lon_sign = E if lon >= 0 else -E
    topos    = wgs84.latlon(lat * N, abs(lon) * lon_sign)
    observer = eph["earth"] + topos

    # Same idea: define widths once, build header and rows with identical spec.
    NW = dict(body=11, date=12, time=5, sign=13, deg=6, type_=14)

    def fmt_nt(body, date, time, sign, deg, type_=''):
        base = (f"  {body:<{NW['body']}} {date:<{NW['date']}} {time:>{NW['time']}}"
                f"  {sign:<{NW['sign']}} {deg:>{NW['deg']}}")
        if type_:
            return base + f"  {type_:<{NW['type_']}}"
        return base

    nt_hdr = fmt_nt('BODY', 'DATE', 'TIME', 'SIGN', 'DEG')
    width  = 70
    bar    = "-" * width

    print()
    print()
    print("=" * width)
    print(f"  NEXT NIGHT PEAK  (scanning up to 366 days from {start_date})")
    print("=" * width)
    print(nt_hdr)
    print(bar)

    for name in ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]:
        sym = PLANET_SYMBOLS[name]
        if name == "Sun":
            print(fmt_nt('Sun', '---', '---', '---', '---', 'always daytime'))
            continue

        if name in ("Venus", "Mercury"):
            result = find_next_night_peak(eph, ts, name, start_date, observer, topos, tz)
            if result is None:
                print(f"  {name:<{NW['body']}} no night visibility in 366 days")
            else:
                sign        = zodiac_sign(result["ecl_lon"])
                deg_in_sign = result["ecl_lon"] % 30
                date_str    = result["date"].strftime("%a %b %d")
                time_str    = fmt_time(result["peak_utc"], tz)
                print(fmt_nt(name, date_str, time_str, sign, f"{deg_in_sign:.1f}", 'night peak'))
            continue

        result = find_next_night_peak(eph, ts, name, start_date, observer, topos, tz)

        if result is None:
            print(f"  {name:<{NW['body']}} no visibility in 366 days")
            continue

        sign        = zodiac_sign(result["ecl_lon"])
        deg_in_sign = result["ecl_lon"] % 30
        date_str    = result["date"].strftime("%a %b %d")
        time_str    = fmt_time(result["peak_utc"], tz)
        type_str    = "transit   " if result["is_transit"] else "night peak"

        print(fmt_nt(name, date_str, time_str, sign, f"{deg_in_sign:.1f}", type_str))

    print(bar)
    print(f"  transit    = upper meridian crossing (highest possible point)")
    print(f"  night peak = highest altitude during darkness (transit is in daytime)")
    print()
    print()
    print("=" * width)
    print()
    print()


# ─────────────────────────────────────────────
#  CONFIGURATION FILE
# ─────────────────────────────────────────────

CONFIG_PATH = pathlib.Path("astro_magick.cfg")

DEFAULT_CONFIG = """[location]
# Latitude in decimal degrees (positive = North, negative = South)
latitude = 37.3382
# Longitude in decimal degrees (positive = East, negative = West)
longitude = -121.8863
# Location name shown in the report header (for your reference only)
name = San Jose, CA

[time]
# IANA timezone name — see https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
timezone = America/Los_Angeles
"""

def load_config():
    """
    Load settings from astro_magick.cfg in the current directory.
    If the file does not exist, write it with the default values and return them.
    Returns (lat, lon, timezone_str).
    """
    if not CONFIG_PATH.exists():
        CONFIG_PATH.write_text(DEFAULT_CONFIG)
        print(f"Created default config file: {CONFIG_PATH.resolve()}", file=sys.stderr)

    cfg = configparser.ConfigParser()
    cfg.read(CONFIG_PATH)

    try:
        lat = cfg.getfloat("location", "latitude")
        lon = cfg.getfloat("location", "longitude")
        tz_str = cfg.get("time", "timezone").strip()
    except (configparser.NoSectionError, configparser.NoOptionError) as e:
        print(f"Error reading config file ({CONFIG_PATH}): {e}", file=sys.stderr)
        sys.exit(1)

    return lat, lon, tz_str


# ─────────────────────────────────────────────
#  ECLIPSE REPORT
# ─────────────────────────────────────────────

LUNAR_ECLIPSE_NAMES = {0: "Penumbral", 1: "Partial", 2: "Total"}

def find_eclipses(eph, ts, start_date, observer, tz):
    """
    Scan forward from start_date up to the DE421 ephemeris limit (2053-10-09).
    Returns a dict with keys:
        lunar_partial  — first penumbral or partial lunar eclipse
        lunar_total    — first total lunar eclipse
        solar_partial  — first partial solar eclipse visible at location
        solar_total    — first total/annular solar eclipse visible at location
    Each value is a dict: {date, time_utc, kind} or None if not found.
    """
    earth    = eph["earth"]
    sun_body = eph["sun"]
    moon_body= eph["moon"]

    # Hard ceiling — DE421 ends 2053-10-09. Leave a 10-day safety margin.
    DE421_LIMIT_TT = ts.utc(2053, 10, 1).tt

    t0 = ts.utc(start_date.year, start_date.month, start_date.day)
    t1 = ts.tt_jd(DE421_LIMIT_TT)

    if t0.tt >= t1.tt:
        return {"lunar_partial": None, "lunar_total": None,
                "solar_partial": None, "solar_total": None}

    results = {
        "lunar_partial": None,
        "lunar_total":   None,
        "solar_partial": None,
        "solar_total":   None,
    }

    # ── Lunar eclipses ────────────────────────────────────────────────────
    try:
        ecl_times, ecl_types, _ = lunar_eclipses(t0, t1, eph)
    except Exception:
        ecl_times, ecl_types = [], []
    for t_ecl, etype in zip(ecl_times, ecl_types):
        kind_str = LUNAR_ECLIPSE_NAMES.get(int(etype), "Unknown")
        entry = dict(
            date     = t_ecl.utc_datetime().date(),
            time_utc = t_ecl.utc_datetime(),
            kind     = kind_str,
        )
        if int(etype) == 2 and results["lunar_total"] is None:
            results["lunar_total"] = entry
        elif int(etype) in (0, 1) and results["lunar_partial"] is None:
            results["lunar_partial"] = entry
        if results["lunar_partial"] and results["lunar_total"]:
            break

    # ── Solar eclipses ────────────────────────────────────────────────────
    def near_new_moon(t):
        try:
            e = earth.at(t)
            sep = e.observe(moon_body).apparent().separation_from(
                  e.observe(sun_body).apparent()).degrees
            return sep < 2.0
        except Exception:
            return False
    near_new_moon.step_days = 25.0

    try:
        events, flags = almanac.find_discrete(t0, t1, near_new_moon)
    except Exception:
        events, flags = [], []

    def sep_at(t):
        obs = observer.at(t)
        return obs.observe(moon_body).apparent().separation_from(
               obs.observe(sun_body).apparent()).degrees

    for t_evt, flag in zip(events, flags):
        if not flag:
            continue
        # Clamp ternary search window strictly within DE421 limit
        lo_tt = max(t0.tt, t_evt.tt - 2)
        hi_tt = min(DE421_LIMIT_TT, t_evt.tt + 2)
        if lo_tt >= hi_tt:
            continue
        try:
            for _ in range(50):
                span = hi_tt - lo_tt
                if span < 1e-6:
                    break
                m1 = ts.tt_jd(lo_tt + span / 3)
                m2 = ts.tt_jd(hi_tt - span / 3)
                if sep_at(m1) > sep_at(m2):
                    lo_tt = lo_tt + span / 3
                else:
                    hi_tt = hi_tt - span / 3
            t_min   = ts.tt_jd((lo_tt + hi_tt) / 2)
            min_sep = sep_at(t_min)
        except Exception:
            continue

        if min_sep > 1.5:
            continue

        is_total = min_sep < 0.3
        kind_str = "Total/Annular" if is_total else "Partial"
        entry = dict(
            date     = t_min.utc_datetime().date(),
            time_utc = t_min.utc_datetime(),
            kind     = kind_str,
        )
        if is_total and results["solar_total"] is None:
            results["solar_total"] = entry
        elif not is_total and results["solar_partial"] is None:
            results["solar_partial"] = entry

        if results["solar_partial"] and results["solar_total"]:
            break

    return results


def print_eclipse_report(eph, ts, start_date, lat, lon, tz):
    """Print next partial and total eclipses (lunar and solar)."""
    lon_sign = E if lon >= 0 else -E
    topos    = wgs84.latlon(lat * N, abs(lon) * lon_sign)
    observer = eph["earth"] + topos

    width = 70
    bar   = "-" * width

    print()
    print()
    print("=" * width)
    print(f"  ECLIPSE FORECAST  (from {start_date} through 2053-10-09)")
    print("=" * width)
    print(f"  {'TYPE':<22} {'DATE':<14} {'TIME (local)':>12}  KIND")
    print(bar)

    eclipses = find_eclipses(eph, ts, start_date, observer, tz)

    rows = [
        ("Lunar  partial",  eclipses["lunar_partial"]),
        ("Lunar  total",    eclipses["lunar_total"]),
        ("Solar  partial",  eclipses["solar_partial"]),
        ("Solar  total",    eclipses["solar_total"]),
    ]
    for label, entry in rows:
        if entry is None:
            print(f"  {label:<22} {'none found before 2053'}")
        else:
            date_str = entry["date"].strftime("%a %b %d %Y")
            time_str = fmt_time(entry["time_utc"], tz)
            print(f"  {label:<22} {date_str:<14} {time_str:>12}  {entry['kind']}")

    print(bar)
    print("  Lunar eclipse times = peak eclipse (geocentric; visible if Moon is up)")
    print("  Solar eclipse times = moment of closest approach at your location")
    print()
    print()
    print("=" * width)
    print()
    print()


# ─────────────────────────────────────────────
#  JSON OUTPUT
# ─────────────────────────────────────────────

def _fmt_local(dt_utc, tz):
    """ISO-8601 local datetime string, or None."""
    if dt_utc is None:
        return None
    return dt_utc.astimezone(tz).strftime("%Y-%m-%dT%H:%M")


def build_daily_json(eph, ts, date, lat, lon, tz):
    """Return a dict with the same data as print_report."""
    lon_sign = E if lon >= 0 else -E
    topos    = wgs84.latlon(lat * N, abs(lon) * lon_sign)
    observer = eph["earth"] + topos

    local_noon = datetime.datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tz)
    t_noon     = ts.from_datetime(local_noon.astimezone(datetime.timezone.utc))

    local_start = datetime.datetime(date.year, date.month, date.day,  0, 0, 0, tzinfo=tz)
    local_end   = datetime.datetime(date.year, date.month, date.day, 23,59,59, tzinfo=tz)
    t0 = ts.from_datetime(local_start.astimezone(datetime.timezone.utc))
    t1 = ts.from_datetime(local_end.astimezone(datetime.timezone.utc))

    sun_body = eph[SKYFIELD_NAMES["Sun"]]
    sun_rise_times, sun_rise_flags = almanac.find_risings(observer,  sun_body, t0, t1, horizon_degrees=-0.8333)
    sun_set_times,  sun_set_flags  = almanac.find_settings(observer, sun_body, t0, t1, horizon_degrees=-0.8333)
    real_sunrises = [t.utc_datetime() for t, f in zip(sun_rise_times, sun_rise_flags) if f]
    real_sunsets  = [t.utc_datetime() for t, f in zip(sun_set_times,  sun_set_flags)  if f]
    sunrise_utc = real_sunrises[0] if real_sunrises else None
    sunset_utc  = real_sunsets[0]  if real_sunsets  else None

    night_start, night_end = night_window(
        sunrise_utc, sunset_utc, date, tz, ts, observer, sun_body)

    def sun_up_at(dt_utc):
        if dt_utc is None:
            return None
        if sunrise_utc is not None and sunset_utc is not None:
            return sunrise_utc <= dt_utc <= sunset_utc
        t_check = ts.from_datetime(dt_utc.replace(tzinfo=datetime.timezone.utc))
        alt, _, _ = observer.at(t_check).observe(sun_body).apparent().altaz()
        return alt.degrees > -0.8333

    phase, illum = get_moon_phase(eph, ts, t_noon)
    eclipse_warnings = check_eclipses(eph, ts, date)

    planets = []
    bodies_lon = {}
    bodies_rx  = {}
    for name in ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]:
        ecl_lon, ra, dec = get_body_position(eph, ts, name, t_noon, observer)
        sign        = zodiac_sign(ecl_lon)
        retrograde  = is_retrograde(eph, ts, name, t_noon)
        deg_in_sign = ecl_lon % 30

        rise_utc, culmination_utc, set_utc, nadir_utc, status = find_events(
            eph, ts, name, date, observer, topos, tz)

        vis_start, vis_end = planet_night_visibility(
            rise_utc, set_utc, status, night_start, night_end)

        planets.append({
            "name":             name,
            "ecliptic_lon":     round(ecl_lon, 4),
            "sign":             sign,
            "deg_in_sign":      round(deg_in_sign, 2),
            "retrograde":       retrograde,
            "status":           status,
            "rise":             _fmt_local(rise_utc, tz),
            "transit":          _fmt_local(culmination_utc, tz),
            "transit_daytime":  sun_up_at(culmination_utc),
            "set":              _fmt_local(set_utc, tz),
            "antitransit":      _fmt_local(nadir_utc, tz),
            "night_vis_start":  _fmt_local(vis_start, tz),
            "night_vis_end":    _fmt_local(vis_end, tz),
        })
        bodies_lon[name] = ecl_lon
        bodies_rx[name]  = retrograde

    conjunctions = [
        {"body_a": a, "body_b": b, "separation_deg": diff,
         "tight": diff <= 3.0}
        for a, b, diff in check_conjunctions(bodies_lon)
    ]

    return {
        "date":              date.isoformat(),
        "location":          {"lat": lat, "lon": lon, "timezone": tz.key},
        "moon_phase":        {"phase": phase, "illumination_pct": illum},
        "eclipse_warnings":  eclipse_warnings,
        "planets":           planets,
        "conjunctions":      conjunctions,
    }


def build_nighttransit_json(eph, ts, start_date, lat, lon, tz):
    """Return a dict with the same data as print_night_transit_report."""
    lon_sign = E if lon >= 0 else -E
    topos    = wgs84.latlon(lat * N, abs(lon) * lon_sign)
    observer = eph["earth"] + topos

    rows = []
    for name in ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]:
        if name == "Sun":
            rows.append({"name": "Sun", "note": "always daytime transit",
                         "date": None, "peak_local": None, "sign": None,
                         "deg_in_sign": None, "type": None})
            continue

        result = find_next_night_peak(eph, ts, name, start_date, observer, topos, tz)
        if result is None:
            rows.append({"name": name, "note": "no visibility found",
                         "date": None, "peak_local": None, "sign": None,
                         "deg_in_sign": None, "type": None})
            continue

        sign        = zodiac_sign(result["ecl_lon"])
        deg_in_sign = round(result["ecl_lon"] % 30, 2)
        peak_local  = _fmt_local(result["peak_utc"], tz)
        type_str    = "transit" if result.get("is_transit") else "night peak"
        note        = "full moon" if result.get("full_moon") else None

        rows.append({
            "name":        name,
            "date":        result["date"].isoformat(),
            "peak_local":  peak_local,
            "sign":        sign,
            "deg_in_sign": deg_in_sign,
            "type":        type_str,
            "note":        note,
        })

    return {
        "start_date": start_date.isoformat(),
        "location":   {"lat": lat, "lon": lon, "timezone": tz.key},
        "planets":    rows,
    }


def build_eclipses_json(eph, ts, start_date, lat, lon, tz):
    """Return a dict with the same data as print_eclipse_report."""
    lon_sign = E if lon >= 0 else -E
    topos    = wgs84.latlon(lat * N, abs(lon) * lon_sign)
    observer = eph["earth"] + topos

    eclipses = find_eclipses(eph, ts, start_date, observer, tz)

    def fmt_entry(e):
        if e is None:
            return None
        return {
            "date":        e["date"].isoformat(),
            "time_local":  _fmt_local(e["time_utc"], tz),
            "time_utc":    e["time_utc"].strftime("%Y-%m-%dT%H:%M") + "Z",
            "kind":        e["kind"],
        }

    return {
        "start_date":     start_date.isoformat(),
        "ephemeris_limit": "2053-10-01",
        "location":       {"lat": lat, "lon": lon, "timezone": tz.key},
        "lunar_partial":  fmt_entry(eclipses["lunar_partial"]),
        "lunar_total":    fmt_entry(eclipses["lunar_total"]),
        "solar_partial":  fmt_entry(eclipses["solar_partial"]),
        "solar_total":    fmt_entry(eclipses["solar_total"]),
    }


def output_json(args, eph, ts, start_date, lat, lon, tz):
    """Collect all requested reports into one JSON object and print it."""
    out = {}

    # Always include daily almanac (one entry per --days)
    daily = []
    for i in range(args.days):
        d = start_date + datetime.timedelta(days=i)
        daily.append(build_daily_json(eph, ts, d, lat, lon, tz))
    out["almanac"] = daily

    if args.nighttransit:
        out["night_transit"] = build_nighttransit_json(eph, ts, start_date, lat, lon, tz)

    if args.eclipses:
        out["eclipses"] = build_eclipses_json(eph, ts, start_date, lat, lon, tz)

    print(json.dumps(out, indent=2, default=str))


# ─────────────────────────────────────────────
#  ENTRY POINT
# ─────────────────────────────────────────────

def main():
    # Load config file first; CLI args override config values when provided.
    cfg_lat, cfg_lon, cfg_tz_str = load_config()

    parser = argparse.ArgumentParser(
        description="✦ Astro Magick — Astronomical almanac for magickal practice",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 astro_magick.py
  python3 astro_magick.py --date 2026-06-21
  python3 astro_magick.py --lat 51.5074 --lon -0.1278 --timezone Europe/London
  python3 astro_magick.py --days 7
  python3 astro_magick.py --lat 89.9 --lon 0.0 --timezone UTC

Config file: astro_magick.cfg (created automatically on first run)
Ephemeris:   de421.bsp (~17MB, downloaded automatically on first run)
        """
    )
    parser.add_argument("--lat",      type=float, default=None,
                        help=f"Latitude  (config: {cfg_lat})")
    parser.add_argument("--lon",      type=float, default=None,
                        help=f"Longitude (config: {cfg_lon})")
    parser.add_argument("--timezone", type=str,   default=None,
                        help=f"IANA timezone name (config: {cfg_tz_str})")
    parser.add_argument("--date",     type=str,   default=None,
                        help="Date YYYY-MM-DD (default: today)")
    parser.add_argument("--days",        type=int,   default=1,
                        help="Number of consecutive days to show (default: 1)")
    parser.add_argument("--nighttransit", action="store_true", default=False,
                        help="Show next night transit for each planet")
    parser.add_argument("--eclipses",     action="store_true", default=False,
                        help="Show next partial and total lunar/solar eclipses")
    parser.add_argument("--json",         action="store_true", default=False,
                        help="Output all requested reports as a single JSON object")
    args = parser.parse_args()

    # CLI args take priority over config file
    lat    = args.lat      if args.lat      is not None else cfg_lat
    lon    = args.lon      if args.lon      is not None else cfg_lon
    tz_str = args.timezone if args.timezone is not None else cfg_tz_str

    try:
        tz = ZoneInfo(tz_str)
    except Exception:
        print(f"Error: unknown timezone '{tz_str}'. Use IANA names like 'America/Los_Angeles'.", file=sys.stderr)
        sys.exit(1)

    if args.date:
        try:
            start_date = datetime.date.fromisoformat(args.date)
        except ValueError:
            print("Error: date must be YYYY-MM-DD format.", file=sys.stderr)
            sys.exit(1)
    else:
        start_date = datetime.date.today()

    print("Loading ephemeris (downloads de421.bsp on first run, ~17MB)...", file=sys.stderr)
    eph, ts = load_ephemeris()

    if args.json:
        output_json(args, eph, ts, start_date, lat, lon, tz)
        return

    for i in range(args.days):
        d = start_date + datetime.timedelta(days=i)
        print_report(eph, ts, d, lat, lon, tz)

    if args.nighttransit:
        print_night_transit_report(eph, ts, start_date, lat, lon, tz)

    if args.eclipses:
        print_eclipse_report(eph, ts, start_date, lat, lon, tz)

if __name__ == "__main__":
    main()
