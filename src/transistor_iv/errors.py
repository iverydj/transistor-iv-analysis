"""Exceptions raised by transistor_iv."""


class CurveError(ValueError):
    """The curve as a whole cannot be analysed (too few points, not a sweep, ...)."""


class ExtractionError(ValueError):
    """One quantity cannot be extracted from this curve (its criterion is never met, ...).

    :func:`transistor_iv.analyze_transfer` catches it per quantity, stores NaN and
    records the message in ``TransferResult.notes``.
    """
