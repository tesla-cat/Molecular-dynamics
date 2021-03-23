
def actionPotential(
  t = None, plot = False,
  period = 4, # ms
):
  offset = -70 # mV
  amp1 = 110 
  if t is None: t = np.linspace(0, period, 100)
  sin = amp1 * np.sin(2*np.pi*(1/period)*t)
  sin[50:] *= 0.5
  sin += offset
  actionPotential = sin
  # E field
  membraneThinkness = 10 # nm
  E = actionPotential / membraneThinkness # mV / nm
  E = 1e-3 / 1e-9 * E # V / m
  if plot:
    plt.plot(t, E)
    plt.xlabel('Time (ms)')
    plt.ylabel('E field (V / m)')
    plt.show()
  return E