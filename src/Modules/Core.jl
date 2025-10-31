module Core

export NumericalWindow, MonoBeam, LaserSource

struct NumericalWindow
  xv
  yv
end

struct LaserSource
  frequency
  power
end

struct MonoBeam
  amplitude
  frequency
end

end
