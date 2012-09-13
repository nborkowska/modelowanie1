#! /usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np

from optparse import OptionParser

class Atom(object):
    
    def __init__(self, position, mass=1.0, energy = 0.0):
        self.position = position
        self.mass = mass
        self.energy = energy                     # potential energy
        self.force = np.zeros(Atoms.dim)
        self.velocity = np.zeros(Atoms.dim)
    
    def getKinEnergy(self):
        return 0.5*self.mass*np.linalg.norm(self.velocity)**2


class RandomSample(object):

    @staticmethod 
    def getSample(mu, sigma, dim):                  
        return 1.5*np.random.normal(mu, sigma, dim)
    
    @staticmethod
    def getSampleSet(mu, sigma, size, dim):
        return [RandomSample.getSample(mu, sigma, dim) for i in xrange(size)]


class Atoms(object):
    """ Container for atoms """
    
    dim = 1                    

    def __init__(self, size):
        self.size = size
        """self.atoms = [Atom(vector) for vector in \
                RandomSample.getSampleSet(0, 10, self.size, self.dim)]"""
        self.atoms = [Atom(np.array([i])) for i in xrange(-self.size/2, self.size/2)]
        #self.atoms = [Atom(np.array([i/10.0])) for i in xrange(self.size)]

    def resetFAndE(self):
        for atom in self.atoms:
            atom.force = np.zeros(self.dim)
            atom.energy = 0


class ForceField(object):
    """ 'Abstract' class"""

    def singleForce(self, atom):
        raise NotImplementedError
    
    def singleEnergy(self, atom):
        raise NotImplementedError
    
    def pairForce(self, atom1, atom2):
        raise NotImplementedError
    
    def pairEnergy(self, atom1, atom2):
        raise NotImplementedError


class SoftWalls(ForceField):
    
    f, L = 0.2, 5                # 10, 10 doesn't work

    def singleForce(self, atom):         
        distance = np.linalg.norm(atom.position)
        if distance < self.L:
            atom.force = np.zeros(Atoms.dim) 
        else:
            atom.force = self.f*(self.L-distance)*atom.position/distance
            

    def singleEnergy(self, atom):
        distance = np.linalg.norm(atom.position)
        if distance < self.L:
            atom.energy = 0.0
        else:
            atom.energy = 0.5*self.f*(self.L - distance)**2
    
    def pairForce(self, atom1, atom2):
        pass
    
    def pairEnergy(self, atom1, atom2):
        pass


class MBM(ForceField):
    
    a, b, c, d = 5.0, 10.0, 3.0, 0.02
    
    def singleForce(self, atom):
        x = atom.position
        atom.force = -self.a*math.exp(-self.b*(x-1)**2)*2*self.b*(x-1) \
                - self.c*math.exp(-(x+1)**2)*2*(1+x) \
                - (x**3)*4*self.d
    
    def singleEnergy(self, atom):
        x = atom.position[0]
        atom.energy = -self.a*math.exp(-self.b*(x-1)**2) \
                - self.c*math.exp(-(x+1)**2) + self.d*x**4

    def pairForce(self, atom1, atom2):
        pass
    
    def pairEnergy(self, atom1, atom2):
        pass


class LenardJones(ForceField):
    
    R, e = 0.5, 1.0

    def singleForce(self, atom):
        pass
    
    def singleEnergy(self, atom):
        pass
    
    def setParams(self, atom1, atom2):
        direction = atom1.position-atom2.position
        distance = np.linalg.norm(atom1.position-atom2.position) # euclidean distance
        a = (self.R/distance)**6
        return direction/distance, a

    def pairForce(self, atom1, atom2):
        normalized, a = self.setParams(atom1, atom2)
        F = -12.0*self.e*a*(a-1)*normalized
        atom1.force += F
        atom2.force -= F

    def pairEnergy(self, atom1, atom2):
        normalized, a = self.setParams(atom1, atom2)
        E = self.e*a*(a-2)
        atom1.energy += E/2
        atom2.energy += E/2


class Simulation(object):
    
    def __init__(self, **kwargs):
        kwargs.update({'system': Atoms(kwargs.pop('noMolecules'))})
        self.__dict__.update(kwargs)

    def trajAndEnergies(self, potential, verlet, previous):
        """ Determine energies and trajectories """
        system = self.system
        trajectory = open('trajectory.xyz', 'w') 
        energies = []
        
        for step in xrange(self.noSteps):
            system.resetFAndE()
            trajectory.write('%d\nkomentarz\n' % (system.size))
            totalPotEnergy = totalKinEnergy = 0.0
            
            for i in xrange(system.size):
                potential.singleForce(system.atoms[i])
                potential.singleEnergy(system.atoms[i])

                for atom in system.atoms[i+1:]:
                    potential.pairForce(system.atoms[i], atom)
                    potential.pairEnergy(system.atoms[i], atom)
                
                totalPotEnergy += system.atoms[i].energy
                totalKinEnergy += system.atoms[i].getKinEnergy()
                trajectory.write('%d\t%d\t0.000\t0.000\n' % (i, system.atoms[i].position[0]))
                previous[i] = verlet.step(system.atoms[i], previous[i], self.stepSize)
            energies.append((totalPotEnergy,totalKinEnergy))

        trajectory.close()
        return energies
    
    def start(self):
        
        prevVel = prevFor = np.zeros((self.system.size,Atoms.dim))
        prev = zip([x.position for x in self.system.atoms], prevVel, prevFor)
        verlet = {'bv': BasicVerlet(), 'vv': VelocityVerlet(), 'lf': LeapFrog()}.get(self.integration)
        potential = {'sw': SoftWalls(), 'mbm': MBM(), 'lj': LenardJones()}.get(self.potential)
        
        means, total = [], []
        energies = self.trajAndEnergies(potential, verlet, prev)
        energy = open('energy.csv', 'w')
        
        """ zapozyczone """
        means.append(np.array(energies[0]))
        energy.write('%f\t%f\t%f\n' % (tuple(means[0])+(sum(means[0]),))) 
        for i in xrange(1,len(energies)):
            means.append(means[i-1]+energies[i])
        for i in xrange(1,len(energies)):
            means[i] /= i+1
            total.append(sum(means[i]))
            energy.write('%f\t%f\t%f\n' % (tuple(means[i])+(sum(means[i]),))) 

        plt.plot(total)
        plt.savefig("overall.svg")
        plt.close()
        energy.close()


class Integration(object):
    """ Abstract """
    
    def step(self, atom, previous, stepSize):
        """ next step for single atom """
        raise NotImplementedError


class BasicVerlet(Integration):
    
    def step(self, atom, previous, stepSize):
        current = (atom.position, atom.velocity, atom.force)
        atom.position = 2*atom.position-previous[0]+(stepSize**2)*atom.force/atom.mass 
        atom.velocity = (atom.position - current[0])/stepSize
        return current
    
    
class VelocityVerlet(Integration):
    
    def step(self, atom, previous, stepSize):
        current = (atom.position, atom.velocity, atom.force)
        atom.velocity = atom.velocity+stepSize*0.5*(previous[2]/atom.mass + atom.force/atom.mass)
        atom.position = atom.position+atom.velocity*stepSize+(stepSize**2)*0.5*atom.force/atom.mass
        return current


class LeapFrog(Integration):
    
    def step(self, atom, previous, stepSize):
        current = (atom.position, atom.velocity, atom.force)
        atom.velocity = atom.velocity+stepSize*atom.force/atom.mass
        atom.position = atom.position+atom.velocity*stepSize
        return current
    

def valid(options):
    """ Check mandatory params """

    for val in options.__dict__.values():
        if val is None:
            return False

    return True

def main(*args):
    
    pOpts = ('sw','mbm','lj')
    iOpts = ('bv','vv','lf')
    
    pHelpText = '%s - soft walls, %s - mbm, %s - Lenard-Jones' % pOpts 
    iHelpText = '%s - basic Verlet, %s - Velocity Verlet, %s - Leapfrog' % iOpts
    
    parser = OptionParser(usage="usage: %prog [options]", version="%prog 1.0")
    parser.add_option('-p', '--potential', type='choice', dest='potential', \
            choices=list(pOpts), help='Used potential, where ' + pHelpText)
    parser.add_option('-i', '--integration', type='choice', dest='integration',
            choices=list(iOpts), help='Integration function, where ' + iHelpText)
    parser.add_option('-n', '--nomolecules', type='int', dest='noMolecules', \
            help='Number of molecules')
    parser.add_option('-s', '--steps', type='int', dest='noSteps', \
            help='Number of simulation steps')
    parser.add_option('-d', '--delta', type='float', dest='stepSize', \
            help='Step length')

    (options, args) = parser.parse_args()

    if not valid(options):
        parser.print_help()
        return 0
       
    if np.__version__ < 1.6:
        print 'Sorry, your numpy version is too old. Please consider upgrading to 1.6.1'
        
    simulation = Simulation(**options.__dict__)
    simulation.start()
    
                                                
if __name__ == '__main__':
    main()
