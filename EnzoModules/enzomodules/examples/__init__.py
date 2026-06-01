"""Worked examples built on the certified EnzoModules kernels.

- ``riemann``: an independent exact Riemann solver (Toro Ch.4), used to
  validate the wrapped kernels against analytic truth.
- ``ppm_sod``: a complete 1D PPM hydro driver that evolves a shock tube using
  only ``enzomodules.hydro.ppm_sweep_1d`` -- the reference for how a downstream
  rewrite would orchestrate the wrapped kernels.
"""
from . import ppm_sod, riemann

__all__ = ["riemann", "ppm_sod"]
