# -*- coding: utf-8 -*-
import re


def parse_slice(spec, multiplier=1):
    """Parses a slice specification like "x" or x:y" or "x:y:step" or "x~y" or "x~y:step".
    The x:y form is an exclusive range, x~y is inclusive.
    Returns a slice object.
    """
    # null spec -- return "all" slice
    if not spec:
        return slice(0, None)
    match = re.match("^(?P<start>\d*)((?P<sep>[:~])(?P<end>\d*)(:(?P<step>\d*))?)?$", spec)
    if not match:
        raise ValueError("invalid slice specification %s" % spec)
    start, end, step = [((match.group(name) or None) and int(match.group(name))) for name in ('start', 'end', 'step')]
    sep = match.group('sep')
    # No match for first separator implies single number. If even that is missing, return full slice
    if sep == None:
        return slice(start * multiplier, (start + 1) * multiplier) if start is not None else slice(0, None)
    # increment end if specified inclusive slice as "start~end"
    elif sep == "~" and end is not None:
        end += 1
    return slice(start and start * multiplier, end and end * multiplier, step and step * multiplier)
