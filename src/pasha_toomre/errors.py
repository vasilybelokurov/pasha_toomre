"""Typed errors raised by the simulation pipeline."""


class PashaToomreError(RuntimeError):
    """Base class for expected project errors."""


class ConfigurationError(PashaToomreError):
    """The requested simulation configuration is invalid."""


class UnknownSofteningError(PashaToomreError):
    """The requested gravitational-softening law is unavailable."""


class DynamicsError(PashaToomreError):
    """The dynamical state is invalid or singular."""


class IntegrationError(PashaToomreError):
    """The numerical integration failed or missed a required event."""


class DiagnosticError(PashaToomreError):
    """A derived orbital diagnostic is not physically or numerically valid."""
