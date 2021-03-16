"""
Use neutron rays simulated with McStas and stored using MCPL_output to simulate
the expected intensity from a GISANS experiment.

For all rays a simulation is conducted and combined in a single intensity graph
binned in wavelength to get actual counts. Resolution is automatically included
as the direction of the events is taken into account.
"""

import mcpl
from numpy import *

def read_events(fname):
    f=mcpl.MCPLFile(fname)
    data=hstack([vstack([pari.weight, pari.x, pari.y, pari.z, pari.ekin, pari.ux, pari.uy, pari.uz]) 
                for pari in f.particle_blocks])
    return data

def calc_aiphil(data):
    l=0.00028601435/sqrt(data[4]) # wavelength in Å from energy in MeV
    ai=arctan2(data[7], -data[5]) # incident angle from direction
    phi=arctan2(-data[6], -data[5]) # incident phi direction
    return ai, phi, l
    
"""
Born Again model to use for events.
"""
import bornagain as ba
from bornagain import deg, angstrom, nm
import numpy

def get_sample():
    """
    Returns a sample with a cylindrically shaped mesocrystal on a substrate.
    """
    # defining materials
    m_ambience = ba.HomogeneousMaterial("Air", 0.0, 0.0)
    m_substrate = ba.HomogeneousMaterial("Substrate", 6e-6, 2e-8)
    m_particle = ba.HomogeneousMaterial("Particle", 6e-4, 2e-8)

    # mesocrystal lattice
    lattice_basis_1 = ba.kvector_t(5.0, 0.0, 0.0)
    lattice_basis_2 = ba.kvector_t(0.0, 5.0, 0.0)
    lattice_basis_3 = ba.kvector_t(0.0, 0.0, 5.0)
    lattice = ba.Lattice(lattice_basis_1, lattice_basis_2, lattice_basis_3)

    # spherical particle that forms the base of the mesocrystal
    sphere_ff = ba.FormFactorFullSphere(2*nm)
    sphere = ba.Particle(m_particle, sphere_ff)

    # crystal structure
    crystal = ba.Crystal(sphere, lattice)

    # mesocrystal
    meso_ff = ba.FormFactorCylinder(20 * nm, 50 * nm)
    meso = ba.MesoCrystal(crystal, meso_ff)

    particle_layout = ba.ParticleLayout()
    particle_layout.addParticle(meso)

    air_layer = ba.Layer(m_ambience)
    air_layer.addLayout(particle_layout)
    substrate_layer = ba.Layer(m_substrate)

    multi_layer = ba.MultiLayer()
    multi_layer.addLayer(air_layer)
    multi_layer.addLayer(substrate_layer)
    return multi_layer


def get_simulation(l, ai, phi):
    """
    Returns a GISAXS simulation with beam and detector defined.
    """
    simulation = ba.GISASSimulation()
    simulation.setDetectorParameters(100, 85.0*deg, 95.0*deg,
                                     50, 0.0*deg, 1.0*deg)
    simulation.setBeamParameters(l*angstrom, ai*deg, phi*deg)
    simulation.getOptions().setIncludeSpecular(True)
    return simulation


def run_simulation(l, ai, phi):
    """
    Runs simulation and returns resulting intensity map.
    """
    simulation = get_simulation(l, ai, phi)
    simulation.setSample(get_sample())
    simulation.runSimulation()
    return simulation.result()

def ba_events(data):
    """
    Generate a 3D intensity array for X,Y,Lambda using a set of McStas events with the BA simulation.
    """
    data_2=calc_aiphil(data)
    data_all=vstack([data, array(data_2)])
    ddata=zeros((50,100,int(data_all[-1].max()*10)+1), dtype=float)
    for i, di in enumerate(data_all.T):
        if i%100==0:
            print(i)
        w=di[0]
        ai,phi,l=di[-3:]
        res=run_simulation(l, ai*180/pi, phi*180/pi).array()
        ddata[:,:,int(l*10)]+=w*res
    return ddata