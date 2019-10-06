import meep as mp
import math
import cmath

base_directory = "/data/home/bt16268/simOutput/pics/"

simulation_size = 20
dpml = 3
sxx = simulation_size + dpml
sxy = simulation_size + dpml
sxz = simulation_size + dpml

resolution = 250
cell_radius = 3
thylakoid_thickness = cell_radius / 3
t = 50

ri_cytoplasm = 1.4
ri_media = 1.33
ri_thylakoid = 1.5

complex_f = True
eps_av = True
symmetry_switch = False


class RunSimulation:

    def __init__(self, dir_name, wavelength, thy_absorb, cyt_absorb):
        self.base_directory = base_directory + str(dir_name)

        self.wavelength = wavelength
        self.thy_absorb = thy_absorb
        self.cyt_absorb = cyt_absorb
        self.frequency = 1 / wavelength
        # Calculate wavelengths dependent on RI
        self.wavelength_in_media = wavelength / ri_media
        self.wavelength_in_cytoplasm = wavelength / ri_cytoplasm
        self.wavelength_in_thylakoid = wavelength / ri_thylakoid
        max_freq = self.frequency - 0.01
        min_freq = self.frequency + 0.01
        self.pulse_width = abs(max_freq - min_freq)

        cell = mp.Vector3(sxx, sxy, 0)
        pml_layers = [mp.PML(dpml)]

        thylakoid_material = mp.Medium(index=ri_thylakoid,
                                       D_conductivity=2 * math.pi * self.frequency * (thy_absorb / (ri_thylakoid ** 2)))
        cytoplasm_material = mp.Medium(index=ri_cytoplasm,
                                       D_conductivity=2 * math.pi * self.frequency * (cyt_absorb / (ri_cytoplasm ** 2)))

        thylakoid_region = mp.Sphere(radius=cell_radius,
                                     center=mp.Vector3(0, 0),
                                     material=thylakoid_material)
        cytoplasm_region = mp.Sphere(radius=cell_radius - thylakoid_thickness,
                                     center=mp.Vector3(0, 0),
                                     material=cytoplasm_material)
        geometry = [thylakoid_region, cytoplasm_region]

        # Sources
        kdir = mp.Vector3(1, 0, 0)  # direction of k (length is irrelevant)
        n = ri_media  # refractive index of material containing the source
        k = kdir.unit().scale(2 * math.pi * self.frequency * n)  # k with correct length

        def pw_amp(k, x0):
            def _pw_amp(x):
                return cmath.exp(1j * k.dot(x + x0))
            return _pw_amp
        source = [
            mp.Source(
                mp.ContinuousSource(frequency=self.frequency, fwidth=self.pulse_width),  # along x axis
                component=mp.Ez,
                center=mp.Vector3(-0.5 * simulation_size, 0),  # x, ,y ,z
                size=mp.Vector3(0, sxy, 0),
                amp_func=pw_amp(k, mp.Vector3(x=-0.5 * simulation_size))
            )
            ]

        sim = mp.Simulation(
            cell_size=cell,
            sources=source,
            boundary_layers=pml_layers,
            resolution=resolution,
            geometry=geometry,
            default_material=mp.Medium(index=ri_media),
            force_complex_fields=complex_f,
            eps_averaging=eps_av,
        )

        sim.use_output_directory(self.base_directory)
        sim.run(mp.at_every(0.6, mp.output_png(mp.Ez, "-Zc/data/home/bt16268/simInput/customColour.txt")), until=t)

        # sim.run(mp.after_time(50,
        #                       mp.at_every(0.6, mp.output_png(
        #                           mp.Ez, "-Zc /data/home/bt16268/simInput/customColour.txt"))), until=t)

# /data/home/bt16268/simInput/customColour.txt


chlorophyll_a = RunSimulation("ChloroA - 430nm Std abs", 0.430, 0.006566, 0)  # 430nm (blue) 652nm (red)
#chlorophyll_b = RunSimulation("700nm no abs", 0.700, 0, 0)  # 453nm (blue) 642nm (red)
#carotenoid = RunSimulation("700nm no abs", 0.700, 0, 0)  # 460nm (blue) 550nm (yellow)

