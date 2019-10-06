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
t = 50

ri_cytoplasm = 1.4
ri_media = 1.33
ri_thylakoid = 1.5

complex_f = True
eps_av = True
symmetry_switch = False


class RunSimulation:

    def __init__(self, dir_name, wavelength, cyt_absorb):
        self.base_directory = base_directory + str(dir_name)

        self.wavelength = wavelength
        self.cyt_absorb = cyt_absorb
        self.frequency = 1 / wavelength
        # Calculate wavelengths dependent on RI
        self.wavelength_in_media = wavelength / ri_media
        self.wavelength_in_cytoplasm = wavelength / ri_cytoplasm
        max_freq = self.frequency - 0.01
        min_freq = self.frequency + 0.01
        self.pulse_width = abs(max_freq - min_freq)

        cell = mp.Vector3(sxx, sxy, 0)
        pml_layers = [mp.PML(dpml)]

        cytoplasm_material = mp.Medium(index=ri_thylakoid,
                                       D_conductivity=2 * math.pi * self.frequency * (cyt_absorb / (ri_thylakoid ** 2)))

        cytoplasm_region = mp.Sphere(radius=cell_radius,
                                     center=mp.Vector3(0, 0),
                                     material=cytoplasm_material)
        geometry = [cytoplasm_region]

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


#UVB_singleCyto = RunSimulation("UV-A 350nm SingleLayer Cytoplasm", 0.350, 0)
#ChloroA_singleCyto = RunSimulation("ChloroA - 430nm Single Layer Cytoplasm", 0.430, 0)
#ChloroB_singleCyto = RunSimulation("ChloroB - 640nm Single Layer Cytoplasm", 0.640, 0)

#UVB_singleThyla = RunSimulation("UV-A 350nm SingleLayer Cytoplasm", 0.350, 0.006566)
#ChloroA_singleThyla = RunSimulation("ChloroA - 430nm Single Layer Cytoplasm", 0.430, 0.006566)
ChloroB_SingleThyla = RunSimulation("ChloroB - 640nm Single Layer Thylakoid", 0.640, 0.006566)


'''

https://watermark.silverchair.com/37-3-395.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAk0wggJJBgkqhkiG9w0BBwagggI6MIICNgIBADCCAi8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMotfnDlbSXNtok5C0AgEQgIICAKfPh7-qMmN1MAGulIEq2mBfkHnMuRy9Y4mIi4_eqEgMQBefe7D3LaGJz5f7WUnfTSnysHwelSx-LKb_yE51BO9lM46pQ5GGisD9eNkVCYxOewmm8HHGdKt4Kc9wAPKSy5AlW1IDb5Vqw6VhXcBbtJcg_r7TOM6sHQZulD2-ssYz0fAZUIXAOoaoh8Cllf5srm3R6rgZC_1IfOgfkPaD0MXixsN2y_T_EJWzftrbCRlUxC8DLmHN0dLDmnpxxwjbvikbzjD9OM1YzgRgqgv1bk3Ax4IGn-aOB5gx4LXlej285-Q_LqrIkxjl4BHyMAA5U4gFekDVE2yi8uZi8cH97VF2yLb5YmaY_QkLlny626BUUujTCAqzVvxuVRgMIUA3omzPCJ42Q_4_WDcqZeQwh2n_GAXIU1d2nG7xqno86FCQSdr7Rc0Pa-53KI9W8r3duWYVcbYdv-DbzQ0D9eRFmanbfsS2yDm9i_IID-MefZ9xFJlPss3cOOx6pZOJ50LL6oY__m1xyS4YCyFq9U4ZHGu5UmYTV6p_C2I30j4Bco7kEg2YSuUP8rfc2nqHh4lad9RoYR964NnEfyUWmCzT4pShFdnEsGULe4NOHD6ycQt664TL7gyANMf_k9USVp81y-KDH2l8sulL7sAOLfwiGUbIc_4Ruuh71HZHTzNZK982

Link to where i found UV-A (320 - 380), UV-B (280 - 320) and UV-C (<280) wavelengths
'''







