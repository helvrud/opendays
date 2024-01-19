#~ import subprocess
import espressomd 
from espressomd import reaction_methods
from espressomd import electrostatics
import os, sys, time, random, pprint, tqdm, py3Dmol, shutil, webbrowser
import numpy as np
pprint = pprint.PrettyPrinter(indent = 4).pprint

class base():
    classname = 'base'
    def __init__(self):
        self.time_step = 0.005
        self.samples = {}

        self.Navogadro = 6.022e23 # 1/mol
        self.kT = 1.38064852e-23*300 # J
        self.RT = self.kT * self.Navogadro # J/mol
        # While there are no interactions sigma is an arbitrary number (could be light-year), but it defines the unit of concentration

        self.epsilon = 1.0 # kT --- a measure of energy
        self.sigma = 1.0 # 0.35nm --- a measure of length
        self.unit_of_length = self.sigma_SI = 0.35 # nm
        self.unit = (self.unit_of_length*1e-9)**3*self.Navogadro*1000 # l/mol
        self.punit = self.kT*self.Navogadro/(self.unit/1000) # J/m3 = Pa

        self.tini = time.time()

        self.eq_steps = 10000 # The number of equilibration steps
        self.N_Samples = 100 # The number of samples


    def set_LJ(self):

        if self.sigma:
            lj_eps = self.epsilon
            lj_sig = self.sigma
            #~ lj_cut = 1.12246
            lj_cut = self.sigma*2**(1./6)
            #~ lj_cut = self.sigma_es+1

            P = list(self.TYPES.keys())
            P.remove('Anode-')
            P.remove('Anode0')
            P.remove('Cathode+')
            P.remove('Cathode0')
            pairs = [(P[i], P[j]) for i in range(len(P)) for j in range(i, len(P))]
            #~ print(particles, pairs)
            self.LJ_params = {}
            for pair in pairs:
                ids = [self.TYPES[pair[0]], self.TYPES[pair[1]]]
                self.system.non_bonded_inter[ids[0], ids[1]].lennard_jones.set_params(epsilon = lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
                self.LJ_params[pair] = self.system.non_bonded_inter[ids[0], ids[1]].lennard_jones.get_params()
                # print(self.system.non_bonded_inter[ids[0], ids[1]].lennard_jones.get_params())
                
            pairs_catplus  = [(P[i], 'Cathode+') for i in range(len(P))]
            pairs_anminus = [(P[i], 'Anode-') for i in range(len(P))]
            pairs_cat0  = [(P[i], 'Cathode0') for i in range(len(P))]
            pairs_an0 = [(P[i], 'Anode0') for i in range(len(P))]
            
            pairs = pairs_catplus+pairs_anminus+ pairs_cat0+pairs_an0
            for pair in pairs:
                ids = [self.TYPES[pair[0]], self.TYPES[pair[1]]]
                self.system.non_bonded_inter[ids[0], ids[1]].lennard_jones.set_params(epsilon = lj_eps, sigma=lj_sig*2, cutoff=2*self.sigma*2**(1./6), shift="auto")
                self.LJ_params[pair] = self.system.non_bonded_inter[ids[0], ids[1]].lennard_jones.get_params()
                # print(self.system.non_bonded_inter[ids[0], ids[1]].lennard_jones.get_params())

    def set_thermostat(self):
        if not hasattr(self, 'seed'):
            self.seed = int(time.time())
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed = self.seed)
        self.thermostat_params = self.system.thermostat.get_state()

    def minimize_energy(self):
        self.system.integrator.set_steepest_descent(f_max=0, gamma=0.1, max_displacement=0.1)
        self.system.integrator.run(100)   # maximal number of steps
        self.system.integrator.set_vv()  # to switch back to velocity Verlet

    def integrate(self):
        #print ('integration '+str(self.steps['md'])+' steps')
        self.system.integrator.run(steps=self.steps['md'])
        
    def reaction(self):
        
        #print ('reaction '+str(self.steps['re'])+' steps')
        self.RE.reaction(steps = self.steps['re'])

    def getN(self, keys = []):
        # ~ if keys == []: keys = self.keys['re']
        if keys == []: keys = self.TYPES.keys()

        for key in keys:
            self.N[key] = len(self.system.part.select(type=self.TYPES[key]))
            try:                self.samples[key] = np.append(self.samples[key], self.N[key])
            except KeyError:    self.samples[key] = np.array([self.N[key]])
        return self.N

    def energy(self):
        energy = self.system.analysis.energy()['total'] - self.system.analysis.energy()['kinetic']
        try:
            self.samples['energy'] = np.append(self.samples['energy'], energy)
        except KeyError:
            self.samples['energy'] = np.array([energy])
        return energy

    def time(self):
        time = self.system.time
        try:
            self.samples['time'] = np.append(self.samples['time'], time)
        except KeyError:
            self.samples['time'] = np.array([time])
        return time

    def pressure(self):
        pressure = self.system.analysis.pressure()['total']
        try:
            self.samples['pressure'] = np.append(self.samples['pressure'], pressure)
        except KeyError:
            self.samples['pressure'] = np.array([pressure])
        return pressure

    def coords(self):
        self.part = {
                'pos': self.system.part.all().pos,
                'id': self.system.part.all().id,
                'type': self.system.part.all().type,
                'bonds': self.system.part.all().bonds,
                'q': self.system.part.all().q,}
        bonds = []
        BONDS = self.part['bonds']
        for i in range(len(BONDS)):
            bonds.append([i])
            for bond in BONDS[i]:
                 bonds[-1].append(bond[1])
            self.part['bonds'] = bonds

        try:             self.Samples['coords'] = np.append(self.Samples['coords'], self.part)
        except KeyError: self.Samples['coords'] = np.array([self.part])

        return self.part

    def writemol2(self):
        #coords = np.random.choice(self.Samples['coords'])
        coords = self.coords()
        # For reference of mol2 file format see the link https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html
        strings = []

        strings.append("# Name: combgel\n")
        strings.append("@<TRIPOS>MOLECULE\n")
        strings.append(self.name.replace('/','_')+"\n")
        strings.append(str(len(coords['id']))+' '+str(len(coords['bonds'])+7)+'\n')
        strings.append("SMALL\n")
        strings.append("USER_CHARGES\n")
        strings.append("@<TRIPOS>ATOM\n")
        
        coords_pos = coords['pos'] % self.box_l
        for i in range(len(coords['id'])):
            idx = coords['id'][i]
            q = coords['q'][i]
            # name = self.NAMES[coords['type'][i]][0].upper()
            name = self.NAMES[coords['type'][i]]
            coord = coords_pos[i]
            #print("atom {} radius 1 name {} type {}\n".format(idx, typ, typ))
            #radius = 0.5
            #if typ in [self.TYPES['nodes'],self.TYPES['PA'],self.TYPES['PHA']]: radius = 1
            strings.append("{} {} {:3.3f} {:3.3f} {:3.3f} {} {} {} {:3.3f}\n".format(idx, name, coord[0], coord[1], coord[2], name, idx, name, q))
        strings.append("@<TRIPOS>BOND\n")
        
        i = 0
        for b in coords['bonds']:
            if len(b) > 0:
                for j in b[1:]:
                    if np.linalg.norm(coords_pos[b[0]] - coords_pos[j])<2.1:
                        strings.append("{} {} {} 1\n".format(i, b[0], j))
                        i+=1
            #print(i)
        strings[3] = str(len(coords['id']))+' '+str(i)+'\n'
        fp = open(self.fname+'.mol', mode='w+t')
        for s in strings:
            fp.write(s)
        fp.close()
        self.fnamemol = self.fname+'.mol'
        print ('% structure saved in '+self.fnamemol)
        return self.fnamemol

##################################################################################
class batt(base):
    classname = 'batt'
    def __init__(self):
        self.name = 'batt'
        #~ self.box_l = box_l
        #~ self.volume = box_l**3
        base.__init__(self)
        # Parameters section
        self.val = {
            'Cl':  -1,
            'Na':  +1,
            'Ca':  +2,
            }

        self.sigma = 1.0 # [ sigma ]
        self.lB = 2.0 # Bjerrum length in water [ sigma ]

        try: idx = max(self.TYPES.values())+1
        except AttributeError: idx = 0; self.TYPES = {}

        type_names = ['Cl', 'Na', 'Ca', 'Anode-','Cathode+','Cathode0','Anode0',]
        for type_name in type_names:
            self.TYPES[type_name] = idx; idx+=1

        self.NAMES = {} # the inverse to TYPES dictionary
        for key in self.TYPES.keys():
            val = self.TYPES[key]
            self.NAMES[val] = key


        self.Samples = {}


        self.keys = {} # This is the dictionary for the key values of interest
        self.keys['md'] = ['pressure']
        self.keys['re'] = []
        


        self.steps = {} # The dictionary storing number of steps for each process ie fhen calling integrate, reaction
        self.steps['md'] = 64
        self.steps['re'] = 32

        self.n_samples = {}
        self.n_samples['re'] = 100 # number of samples in one reaction run step
        self.n_samples['md'] = 50 # number of samples in one md step
 
        self.exclusion_radius = 1.0*self.sigma

    def set_electrodes_reaction(self,charge = 0):
        
        if abs(self.pKe) != np.infty: 
            if not hasattr(self,'RE'):
                self.RE = reaction_methods.ReactionEnsemble(
                        kT=1., 
                        exclusion_range=0.5, 
                        seed = np.random.randint(100))
            
            self.Gamma = self.unit*10**(-self.pKe)
            self.RE.add_reaction(gamma = self.Gamma, 
                reactant_types=[self.TYPES['Cathode+'],self.TYPES['Anode-']],
                reactant_coefficients=[1,1],
                product_types=[self.TYPES['Cathode0'],self.TYPES['Anode0']],
                product_coefficients=[1,1],
                default_charges={self.TYPES['Cathode+']:+charge, self.TYPES['Anode-']:-charge, self.TYPES['Cathode0']:0, self.TYPES['Anode0']:0, })

            print ('Setting up electrode reaction , Gamma = '+str(self.Gamma))
            self.RE_params = self.RE.get_status()
            self.RE_params['volume'] = self.RE.get_volume()

    def show(self):
        fnamemol = self.writemol2()
        shutil.copyfile(fnamemol, '/home/kvint/web/'+fnamemol)
        webbrowser.open('127.0.0.1/mol.html?'+fnamemol)
        
        
    def __str__(self):
        self.name = self.classname
        # This generates the suffixes by the dictionary self.p

        self.name += '_box_l'+str(np.round(self.box_l,3))
        if self.lB != 2.0:
            self.name += '_lB'+str(self.lB)
        if self.sigma != 1.0:
            self.name += '_sigma'+str(self.sigma)
        if self.epsilon != 1.0:
            self.name = self.name+'_epsilon'+str(self.epsilon)
        if self.N_Samples != 100: self.name += '_N'+str(self.N_Samples)

        self.fname = 'data/'+self.name
        self.fnameout = self.fname+".out"
        self.fnamepkl = self.fname+".pkl"
        self.fnamepy = self.fname+".py"

        # ~ self.fnamerun     = self.fname+'.run'
        self.fnameqsub    = self.fname+'.qsub'
        self.fnameqsubout = self.fname+'.qsubout'
        self.fnameqsuberr = self.fname+'.qsuberr'

        return self.name

    def update_re_samples(self):
        return self.getN(keys = self.keys['re'])
        
    def update_md_samples(self):
        ''' returns either mindist or pressure(if sigma is not zero)'''
        fl = (not self.sigma)*[self.mindist]+bool(self.sigma)*[self.pressure]
        return [f() for f in fl]

    def set_EL(self,prefactor=1.0, gap_size = 0):
        if self.lB:
            print('\n #### Seting the Electrostatics ####')
            Q = self.system.part.all().q

            tini = time.time();
            if Q.any():
                self.p3m = espressomd.electrostatics.P3M(prefactor=self.lB, accuracy=1e-4)
                if gap_size: 
                    #elc = espressomd.electrostatics.ELC(actor=p3m, gap_size=gap_size, maxPWerror=1e-3)
                    # This simulates the metalic electrodes at z = 0 and z = Lz -h
                    self.elc = espressomd.electrostatics.ELC(actor=self.p3m, gap_size=gap_size, maxPWerror=1e-3, const_pot=True, delta_mid_bot=1)
                    #self.system.actors.add(self.elc)
                    self.p3m_params = self.elc.get_params()
                else:
                    #self.system.actors.add(self.p3m)
                    self.system.electrostatics.solver = self.p3m
                    self.p3m_params = self.p3m.get_params()
                #self.p3m.tune(accuracy=1e-4)
                
                print("Output of p3m captured.")
                uptime = time.time() - tini
                print('Done, uptime = ', uptime)
            else:
                print('No charges --- no P3M, uptime = ')
   

   
################################################################################
b = batt()
cs = 0.1 # mol/l --- salinity of electrolyte
rho = 0.5 # density of charges on the electrode


b.qsu = 1 # d


b.N_Samples = 100 # A number of samples generated by sample() function
b.sigma = 1.0
#s.lB = 2.0
b.lB = 0.0


cna = cs


nna = 200 # A number of sodium ions
nca = int(np.ceil(nna*0.114)) # A number of calcium ions (a fraction of nna)
ncl = nna + 2*nca #number of chloride ions (neutralizing Na and Ca)

cs_es = cs*b.unit # A desired density of the Na ions in a simulation box 
b.box_l = (nna/cs_es)**(1./3) # The zize of the box which guarantee the desired cs_es
b.volume = b.box_l**3 #Volime of the box
print ('box_l = ', b.box_l)

b.system = espressomd.System(box_l = 3*[b.box_l])

type_suplus = b.TYPES['Cathode+']
type_suminus = b.TYPES['Anode-']
type_cl = b.TYPES['Cl']
type_na = b.TYPES['Na']
type_ca = b.TYPES['Ca']


# Add surface
Nsu_ = int(np.ceil(b.box_l/b.sigma)) # number of beads in one dimension
Nsu = Nsu_**2                        # number of beads in the electrode
qsu = int(rho*b.box_l**2) / Nsu      # the valence of a bead of electrode
# qsu = -1
b.gap_size = b.box_l/3

# First surface
for i in range(Nsu_):
    for j in range(Nsu_):
        pos=[i,j,b.gap_size]
        b.system.part.add(pos=pos, type=b.TYPES['Cathode+'], q=qsu, fix = 3*[True])

# Second surface
for i in range(Nsu_):
    for j in range(Nsu_):
        pos=[i,j, 2 * b.gap_size]
        b.system.part.add(pos=pos, type=b.TYPES['Anode-'], q=-qsu, fix = 3*[True])

b.qsu = qsu

for i in range(nna):
    pos = np.random.rand(3)*b.box_l*[1,1,1/3] + [0,0,b.box_l/3]
    b.system.part.add(pos=pos, type=type_na, q=1)

for i in range(ncl):
    pos = np.random.rand(3)*b.box_l*[1,1,1/3] + [0,0,b.box_l/3]
    b.system.part.add(pos=pos, type=type_cl, q=-1)

for i in range(nca):
    pos = np.random.rand(3)*b.box_l*[1,1,1/3] + [0,0,b.box_l/3]
    b.system.part.add(pos=pos, type=type_ca, q=2)

b.system.time_step = b.time_step
b.system.cell_system.skin = 0.4

b.integrate()
b.set_LJ()
print(b)




   
