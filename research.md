---
layout: default
---

# Research
---
Full lists of my publications can be found on [ADS](https://ui.adsabs.harvard.edu/#search/q=author%3A%22Halevi%2C%20Goni%22&sort=date%20desc%2C%20bibcode%20desc)/[ORCID](https://orcid.org/0000-0002-7232-101X)/[arXiv](https://arxiv.org/search/?query=halevi%2Cgoni&searchtype=all&source=header)

## If you are interested in hearing about what I'm working on these days, reach out and invite me to give a talk at your institution!

---
<!---
## Works in progress

### X-ray search for dwarf galaxies hosting accreting black holes
I work with [Jenny Greene](https://crispygreene.wixsite.com/jenny) and [Andy Goulding](https://www.astro.princeton.edu/~goulding/) to find X-ray bright AGN in low-mass (_M_<sub>*</sub> < 10<sup>10</sup> _M_<sub>☉</sub>) galaxies out to _z_ ~ 1. We use Subaru data from the [Hyper Suprime-Cam Survey](https://hsc.mtk.nao.ac.jp/ssp/) with additional coverage by the CFHT through [CLAUDS](http://www.ap.smu.ca/~agolob/clauds/survey/). These data provide us with accurate photometric redshifts, and thus stellar masses, in a wide and deep region of the sky. We crossmatch the positions of low-mass galaxies with XMM Newton data in the XMM-LSS field and perform various techniques designed to confirm the nature of the X-ray sources. In the end, we build up a sample of roughly 80 accreting massive black holes in low-mass hosts. 

<center><img src="figures/M_z_AGN.png" alt="mass-redshift" width="700px"/></center>
<center><small>[Preliminary results] Stellar mass as a function of redshift for our X-ray selected AGN sample. Color indicates hard-band X-ray luminosities.</small></center> 


### Nuclear reactions in NS-WD accretion disks
I have been working to implement a nuclear reaction network within the MHD code [Athena++](https://princetonuniversity.github.io/athena/) with the goal of simulating the evolution of the torus that forms when a white dwarf is disrupted by a companion neutron star or stellar mass black hole. Accretion disks formed in this way are interesting because the timescales for nuclear reactions to occur can be comparable to the dynamical time in the disk, such that energy deposited by these reactions can dynamically affect the disk structure. Past work has shown that these events power a relatively low-luminosity, rapidly-evolving transient. The hydrodynamical 2D simulations of [Fernández, Margalit, & Metzger (2019)](https://ui.adsabs.harvard.edu/abs/2019arXiv190506343F/abstract) show that nuclear reactions cause the disk to exhibit an "onion-like" structure in composition, as predicted by the early 1D, steady-state calculations of [Metzger (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.419..827M/abstract). As part of my thesis work with [Jim Stone](https://www.astro.princeton.edu/~jstone/), I plan to investigate the effects of including magnetic fields in simulations of such systems, and thus having real MRI-driven turbulence rather than an unphysical alpha viscosity. Other relevant references include: [Zenati et al. (2018)](https://ui.adsabs.harvard.edu/#abs/2018arXiv180709777Z/abstract); [Fernández & Metzger (2013)](https://ui.adsabs.harvard.edu/#abs/2013ApJ...763..108F/abstract); and [Margalit & Metzger (2016)](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.461.1154M/abstract).

<center><img src="figures/onion.png" alt="FM2013" width="700px"/></center>
<center><small>Time-averaged onion-like structure in the disk and outflows produced by the disruption of a carbon/oxygen white dwarf by a neutron star. Source: Fernández, Margalit, & Metzger (2019).</small></center>

### Modeling type I X-ray bursts
Since starting at Princeton in September 2017, I have been running simplified simulations of Type I X-ray bursts and modeling their light curves. These bursts occur when an accreting neutron star in a low-mass X-ray binary (LMXB) is ignited due to the compression of accreted material at the NS surface. Observations of these bursts leave some questions unanswered. For example, it is unclear where the ignition first occurs, how fast it spreads, and how these factors affect the observed light curves. There are also mysteries such as why the high frequency oscillations in the light curve, presumably caused by the rotation of the star, get higher in frequency as the burst proceeds. In order to shed light on some of this uncertainty, I execute spectral codes ([SPHEREPACK](https://www2.cisl.ucar.edu/resources/legacy/spherepack) and [Dedalus](http://dedalus-project.org)) that solve the two-dimensional shallow-water equations to evolve the burning front on the spherical surfaces of model neutron stars. I am experimenting with different parameters including the timescale for burning and the latitude of ignition. I also wrote a light curve generator that takes into account general relativistic effects in calculating the resultant light curves from these simulations. The goal is to compare simulated light curves with observed ones and draw some conclusions about the physics of the ignition and burning front propagation in observed bursts. This work, done in collaboration with [Anatoly Spitkovsky](https://www.astro.princeton.edu/~anatoly/), largely builds off of [Spitkovsky et al. (2002)](https://ui.adsabs.harvard.edu/#abs/2002ApJ...566.1018S/abstract).

<center><img src="figures/nsburn.png" alt="nsburn" width="700px"/></center>
<center><small>[Preliminary results] Simulated light curves for the rise of an X-ray burst in the case where the nuclear burning time is 50 ms and the ignition occurs at a latitude of 75 degrees. The different colors represent different viewing inclinations.</small></center>

### Dynamics of misaligned jet-driven supernovae
I continue to work with [Philipp Mösta](https://pmoesta.wordpress.com/) on understanding the dynamics of the simulations we ran for the project described below. In particular, we ran several models in which the initial core had a magnetic field axis that was misaligned from the rotation axis. This misalignment yields interesting MHD dynamics that we are currently exploring.

---

## Completed (undergraduate) work
### r-process nucleosynthesis in jet-driven supernovae
Working in collaboration with [Philipp Mösta](https://pmoesta.wordpress.com/) and others, beginning in early 2016, I explored the question of where the heaviest elements in the universe are synthesized in space. The conditions required for the r-process, a rapid succession of neutron captures followed by beta decays, are rare. The question of where, astrophysically, these conditions are met is one that people have been trying to answer for decades without much success. Recently, and especially in the past year, a lot of progress has been made; in particular, it's been shown directly that r-process elements get produced in neutron star mergers. However, there remain details to be worked out: are they produced anywhere else in the universe? What are the typical yields and rates of the events that produce r-process elements? Are they compatible with observations of abundances in metal-poor stars and dwarf galaxies? It remains possible that there is another channel, or even multiple channels, for r-process nucleosynthesis to occur, although it may not dominate in most environments.

To investigate this, we specifically looked at the possibility of a special class of core-collapse supernovae producing r-process elements. These supernovae require rapidly rotating and highly magnetized pre-collapse cores that produce jetted explosions. Our simulations of these events were run with the general relativistic 3D MHD [Einstein Toolkit](https://einsteintoolkit.org/) and in post-processing, we used a nuclear reaction network called [SkyNet](http://jonaslippuner.com/research/skynet/) to calculate the final abundances of each isotope. We performed a robust study that included various magnetic field orientations, parameterization of neutrino effects, and evidence of the importance of 3D dynamics. Ultimately, we found that with sufficiently strong initial magnetic fields, the thermodynamic conditions in the jet and in particular the abundance of free neutrons are sufficient for producing the full spectrum of r-process elements. However, the magnetic fields that we need to impose are potentially unrealistically strong, calling into question the applicability of this channel in the real universe.
For more information, see the two papers that have come out of these endeavours: [**Halevi** & Mösta (2018)](https://ui.adsabs.harvard.edu/#abs/2018MNRAS.477.2366H/abstract) and [Mösta, Roberts, **Halevi** et al. (2018)](https://ui.adsabs.harvard.edu/#abs/2017arXiv171209370M/abstract).

<center><img src="figures/tracers.png" alt="tracers" height="450px"/> <img src="figures/volren.png" alt="volren" height="450px"/></center>
<center><small>For an initial field of 10<sup>13</sup> G, we get a jetted explosion that produces r-process elements. Here, we show the tracer particle trajectories for a small sample of the particles, and a volume rendering colored by specific entropy.</small></center>

### Observations: supernovae and active galactic nuclei
As an undergrad, I mostly did observational work as part of [Alex Filippenko](http://w.astro.berkeley.edu/~alex/)'s research team. I discovered several supernovae by checking images taken by the group's robotic telescope, the Katzmann Automated Imaging Telescope ([KAIT](http://w.astro.berkeley.edu/bait/kait.html)), which is located on Mount Hamilton at [Lick Observatory](http://mthamilton.ucolick.org/). Additionally, I spent many nights taking photometric data, mostly of recently detected SNe, with the Nickel, a 1-m telescope at Lick. I also spent many nights acquiring spectra of mainly SNe and some cepheids and AGN with Shane, the observatory's 3-m telescope. My AGN spectra were used for the Lick AGN Mapping Project ([LAMP](https://www.physics.uci.edu/~barth/lamp.html)), a collaborative reverberation mapping campaign led by [Aaron Barth](https://www.physics.uci.edu/~barth/).

-->
