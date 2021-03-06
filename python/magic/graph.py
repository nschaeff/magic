# -*- coding: utf-8 -*-
import numpy as np
import os
import re
from .log import MagicSetup
from .libmagic import scanDir
from magic.setup import buildSo
from .npfile import npfile

if buildSo:
    try:
        import sys
        if sys.version_info.major == 3:
            import magic.greader_single3 as Gsngl
            import magic.greader_double3 as Gdble
        elif sys.version_info.major == 2:
            import magic.greader_single2 as Gsngl
            import magic.greader_double2 as Gdble
        readingMode = 'f2py'
    except ImportError:
        readingMode = 'python'
    # print('read with %s' % readingMode)
else:
    readingMode = 'python'


def getGraphEndianness(filename):
    """
    This function determines the endianness of the graphic files

    :param filename: input of the filename
    :type filename: str
    :returns: the endianness of the file ('B'='big_endian'
              or 'l'='little_endian')
    :rtype: str
    """
    f = npfile(filename, endian='B')
    try:
        f.fort_read('S20')
        endian = 'B'
    except TypeError:
        endian = 'l'
    f.close()

    return endian


class MagicGraph(MagicSetup):
    """
    This class allows to read the 3-D graphic outputs of the MagIC code
    (:ref:`G_#.TAG <secGraphFile>` and G_ave.TAG) files. Those are
    binary unformatted outputs, there are therefore two ways to load them:

       * If buildLib=True in magic.cfg and the fortran libraries were correctly
         built, then the reader uses a fortran program that is expected to be
         much faster than the pure python routine.
       * If buildLib=False, then a pure python program is used to read the G
         files.

    >>> # Regular G files
    >>> gr = MagicGraph(ivar=1, tag='N0m2a')
    >>> print(gr.vr.shape) # shape of vr
    >>> print(gr.ek) # print ekman number
    >>> print(gr.minc) # azimuthal symmetry
    >>> # Averaged G file with double precision
    >>> gr = MagicGraph(ave=True, tag='N0m2', precision='Float64')
    """

    def __init__(self, ivar=None, datadir='.', quiet=True,
                 ave=False, tag=None, precision='Float32'):
        """
        :param ave: when set to True, it tries to find an average
                    G file (G_ave.TAG)
        :type ave: bool
        :param ivar: the number of the G file
        :type ivar: int
        :param tag: extension TAG of the G file. If not specified, the
                    most recent G_#.TAG file found in the directory will
                    be selected.
        :type tag: str
        :param quiet: when set to True, makes the output silent
        :type quiet: bool
        :param datadir: directory of the G file (default is . )
        :type datadir: str
        :param precision: single or double precision (default 'Float32')
        :type precision: str
        """
        self.precision = precision

        if ave:
            self.name = 'G_ave'
        else:
            self.name = 'G_'

        if tag is not None:
            if ivar is not None:
                file = '%s%i.%s' % (self.name, ivar, tag)
                filename = os.path.join(datadir, file)
            else:
                files = scanDir('%s*%s' % (self.name, tag))
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists('log.%s' % tag):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ivar is not None:
                files = scanDir('%s%i*' % (self.name, ivar))
                filename = os.path.join(datadir, files[-1])
            else:
                files = scanDir('%s*' % self.name)
                filename = os.path.join(datadir, files[-1])
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists('log.%s' % ending):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        # Get file endianness
        endian = getGraphEndianness(filename)

        if readingMode != 'python':

            if self.precision == 'Float32':
                G = Gsngl.greader_single
            elif self.precision == 'Float64':
                G = Gdble.greader_double

            G.readg(filename, endian=endian)
            self.nr = int(G.nr)
            self.nr_ic = int(G.nric)-1
            self.ntheta = int(G.nt)
            self.npI = int(G.np)
            self.minc = int(G.minc)
            self.time = G.time
            self.ra = G.ra
            self.ek = G.ek
            self.pr = G.pr
            self.prmag = G.prmag
            self.radratio = G.radratio
            self.sigma = G.sigma
            if self.npI == self.ntheta*2:
                self.npI = int(self.npI/self.minc)
            self.nphi = int(self.npI*self.minc+1)
            self.radius = G.radius
            self.colatitude = G.colat
            self.entropy = G.entropy
            self.vr = G.vr
            self.vtheta = G.vt
            self.vphi = G.vp
            self.pre = G.pre
            self.xi = G.xi
            if self.prmag != 0:
                self.Br = G.br
                self.Btheta = G.bt
                self.Bphi = G.bp

            if self.prmag != 0 and self.nr_ic > 1:
                self.radius_ic = G.radius_ic
                self.Br_ic = G.br_ic
                self.Btheta_ic = G.bt_ic
                self.Bphi_ic = G.bp_ic
        else:
            # Read data
            inline = npfile(filename, endian=endian)

            # Read the header
            version = inline.fort_read('|S20')[0]
            version = version.rstrip()
            runID = inline.fort_read('|S64')[0]
            self.time, n_r_max, n_theta_max, self.npI, n_r_ic_max, minc, \
                nThetaBs, self.ra, self.ek, self.pr, self.prmag, \
                self.radratio, self.sigma = \
                inline.fort_read(self.precision)

            self.nr = int(n_r_max)
            self.ntheta = int(n_theta_max)
            self.npI = int(self.npI)
            self.minc = int(minc)
            if self.npI == self.ntheta*2:
                self.npI = int(self.npI/self.minc)
            self.nphi = self.npI*self.minc+1
            self.nr_ic = int(n_r_ic_max) - 1
            self.nThetaBs = int(nThetaBs)

            if not quiet:
                print('Rayleigh = %.1e, Ekman = %.1e, Prandtl = %.1e' %
                      (self.ra, self.ek, self.pr))
                print('nr = %i, nth = %i, nphi = %i' %
                      (self.nr, self.ntheta, self.npI))

            self.colatitude = inline.fort_read(self.precision)
            self.radius = np.zeros((self.nr), self.precision)

            if self.prmag != 0 and self.nr_ic > 1:
                self.radius_ic = np.zeros((self.nr_ic), self.precision)

            entropy = np.zeros((self.npI, self.ntheta, self.nr),
                               self.precision)
            vr = np.zeros_like(entropy)
            vtheta = np.zeros_like(entropy)
            vphi = np.zeros_like(entropy)
            if self.prmag != 0:
                Br = np.zeros_like(entropy)
                Btheta = np.zeros_like(entropy)
                Bphi = np.zeros_like(entropy)
                if self.sigma != 0:
                    Br_ic = np.zeros((self.npI, self.ntheta, self.nr_ic),
                                     self.precision)
                    Btheta_ic = np.zeros_like(Br_ic)
                    Bphi_ic = np.zeros_like(Br_ic)
            if version == b'Graphout_Version_8' or \
               version == b'Graphout_Version_10' or \
               version == b'Graphout_Version_12':
                pressure = np.zeros_like(entropy)
            if version == b'Graphout_Version_11' or \
               version == b'Graphout_Version_12':
                xi = np.zeros_like(entropy)

            for k in range(self.nr*self.nThetaBs):
                # radius and Thetas in this block
                ir, rad, ilat1, ilat2 = inline.fort_read(self.precision)
                ir = int(ir)
                ilat1 = int(ilat1) - 1
                ilat2 = int(ilat2) - 1
                self.radius[ir] = rad
                nth_loc = ilat2 - ilat1 + 1

                if version == b'Graphout_Version_9' or \
                   version == b'Graphout_Version_10' or \
                   version == b'Graphout_Version_11' or \
                   version == b'Graphout_Version_12':
                    if self.prmag != 0:
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        entropy[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        vr[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        vtheta[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        vphi[:, ilat1:ilat2+1, ir] = data.T
                        if version == b'Graphout_Version_11' or \
                           version == b'Graphout_Version_12':
                            data = inline.fort_read(self.precision,
                                                    shape=(nth_loc, self.npI))
                            xi[:, ilat1:ilat2+1, ir] = data.T
                        if version == b'Graphout_Version_10' or \
                           version == b'Graphout_Version_12':
                            data = inline.fort_read(self.precision,
                                                    shape=(nth_loc, self.npI))
                            pressure[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        Br[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        Btheta[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        Bphi[:, ilat1:ilat2+1, ir] = data.T

                    else:
                        # vectorize
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        entropy[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        vr[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        vtheta[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        vphi[:, ilat1:ilat2+1, ir] = data.T
                        if version == b'Graphout_Version_11' or \
                           version == b'Graphout_Version_12':
                            data = inline.fort_read(self.precision,
                                                    shape=(nth_loc, self.npI))
                            xi[:, ilat1:ilat2+1, ir] = data.T
                        if version == b'Graphout_Version_10' or \
                           version == b'Graphout_Version_12':
                            data = inline.fort_read(self.precision,
                                                    shape=(nth_loc, self.npI))
                            pressure[:, ilat1:ilat2+1, ir] = data.T
                else:
                    if self.prmag != 0:
                        # To vectorize, one must also read the end-of-line
                        # symbol that accounts for 2 additionnal floats
                        # For the last line of the chunk, we should not
                        # read the last symbol
                        if version == b'Graphout_Version_8':
                            data = inline.fort_read(
                                self.precision,
                                shape=(8*(self.npI+2)*nth_loc-2))
                            # Add 2 zeros for the last lines
                            data = np.append(data, 0.)
                            data = np.append(data, 0.)
                            data = data.reshape((8, nth_loc, self.npI+2))
                        else:
                            data = inline.fort_read(
                                self.precision,
                                shape=(7*(self.npI+2)*nth_loc-2))
                            # Add 2 zeros for the last lines
                            data = np.append(data, 0.)
                            data = np.append(data, 0.)
                            data = data.reshape((7, nth_loc, self.npI+2))
                        data = data[:, :, :-2:]
                        entropy[:, ilat1:ilat2+1, ir] = data[0, ...].T
                        vr[:, ilat1:ilat2+1, ir] = data[1, ...].T
                        vtheta[:, ilat1:ilat2+1, ir] = data[2, ...].T
                        vphi[:, ilat1:ilat2+1, ir] = data[3, ...].T
                        if version == b'Graphout_Version_8':
                            pressure[:, ilat1:ilat2+1, ir] = data[4, ...].T
                            Br[:, ilat1:ilat2+1, ir] = data[5, ...].T
                            Btheta[:, ilat1:ilat2+1, ir] = data[6, ...].T
                            Bphi[:, ilat1:ilat2+1, ir] = data[7, ...].T
                        else:
                            Br[:, ilat1:ilat2+1, ir] = data[4, ...].T
                            Btheta[:, ilat1:ilat2+1, ir] = data[5, ...].T
                            Bphi[:, ilat1:ilat2+1, ir] = data[6, ...].T
                    else:
                        if version == b'Graphout_Version_8':
                            data = inline.fort_read(
                                self.precision,
                                shape=(5*(self.npI+2)*nth_loc-2))
                            # Add 2 zeros for the last lines
                            data = np.append(data, 0.)
                            data = np.append(data, 0.)
                            data = data.reshape((5, nth_loc, self.npI+2))
                        else:
                            data = inline.fort_read(
                                self.precision,
                                shape=(4*(self.npI+2)*nth_loc-2))
                            # Add 2 zeros for the last lines
                            data = np.append(data, 0.)
                            data = np.append(data, 0.)
                            data = data.reshape((4, nth_loc, self.npI+2))
                        data = data[:, :, :-2:]
                        entropy[:, ilat1:ilat2+1, ir] = data[0, ...].T
                        vr[:, ilat1:ilat2+1, ir] = data[1, ...].T
                        vtheta[:, ilat1:ilat2+1, ir] = data[2, ...].T
                        vphi[:, ilat1:ilat2+1, ir] = data[3, ...].T
                        if version == b'Graphout_Version_8':
                            pressure[:, ilat1:ilat2+1, ir] = data[4, ...].T

            if self.prmag != 0 and self.nr_ic > 1:
                for k in range(self.nr_ic):
                    # radius and Thetas in this block
                    ir, rad, ilat1, ilat2 = inline.fort_read(self.precision)
                    ir = ir-self.nr
                    ir = int(ir)
                    ilat1 = int(ilat1) - 1
                    ilat2 = int(ilat2) - 1
                    self.radius_ic[ir] = rad
                    nth_loc = ilat2 - ilat1 + 1

                    if version == b'Graphout_Version_9' or \
                       version == b'Graphout_Version_10' or \
                       version == b'Graphout_Version_11' or \
                       version == b'Graphout_Version_12':

                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        Br_ic[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        Btheta_ic[:, ilat1:ilat2+1, ir] = data.T
                        data = inline.fort_read(self.precision,
                                                shape=(nth_loc, self.npI))
                        Bphi_ic[:, ilat1:ilat2+1, ir] = data.T

                    else:
                        data = inline.fort_read(self.precision,
                            shape=(3*(self.npI+2)*nth_loc-2))
                        data = np.append(data, 0.)
                        data = np.append(data, 0.)
                        data = data.reshape((3, nth_loc, self.npI+2))
                        data = data[:, :, :-2:]
                        Br_ic[:, ilat1:ilat2+1, ir] = data[0, ...].T
                        Btheta_ic[:, ilat1:ilat2+1, ir] = data[1, ...].T
                        Bphi_ic[:, ilat1:ilat2+1, ir] = data[2, ...].T

            inline.close()
            # Sorting of data (strange hemispherical way that the magic
            # code uses
            entropy = self.rearangeLat(entropy)
            vr = self.rearangeLat(vr)
            vtheta = self.rearangeLat(vtheta)
            vphi = self.rearangeLat(vphi)
            if self.prmag != 0:
                Br = self.rearangeLat(Br)
                Btheta = self.rearangeLat(Btheta)
                Bphi = self.rearangeLat(Bphi)
                if self.nr_ic > 1:
                    Br_ic = self.rearangeLat(Br_ic)
                    Btheta_ic = self.rearangeLat(Btheta_ic)
                    Bphi_ic = self.rearangeLat(Bphi_ic)
            if version == b'Graphout_Version_11' or \
               version == b'Graphout_Version_12':
                xi = self.rearangeLat(xi)
            if version == b'Graphout_Version_8' or \
               version == b'Graphout_Version_10' or \
               version == b'Graphout_Version_12':
                pressure = self.rearangeLat(pressure)

            # Normalise r
            self.radius = self.radius/(1.-self.radratio)
            if self.prmag != 0 and self.nr_ic > 1:
                self.radius_ic = self.radius_ic/(1.-self.radratio)

            # Full solution by repeating minc times the structure
            self.nphi = self.minc * self.npI + 1
            self.entropy = entropy
            self.vr = vr
            self.vtheta = vtheta
            self.vphi = vphi
            if self.prmag != 0:
                self.Br = Br
                self.Btheta = Btheta
                self.Bphi = Bphi
                if self.nr_ic > 1:
                    self.Br_ic = Br_ic
                    self.Btheta_ic = Btheta_ic
                    self.Bphi_ic = Bphi_ic

            if version == b'Graphout_Version_11' or \
               version == b'Graphout_Version_12':
                self.xi = xi
            if version == b'Graphout_Version_8' or \
               version == b'Graphout_Version_10' or \
               version == b'Graphout_Version_12':
                self.pre = pressure

    def rearangeLat(self, field):
        """
        This function is used to unfold the colatitudes

        :param field: input array with MagIC ordering of colatitudes (i.e.
                      successively Northern Hemisphere and Southern
                      Hemisphere)
        :type field: numpy.ndarray
        :return: an array with the regular ordering of the colatitudes
        :rtype: numpy.ndarray
        """
        even = field[:, ::2, :]
        odd = field[:, 1::2, :]

        return np.concatenate((even, odd[:, ::-1, :]), axis=1)
