#!/usr/bin/env python
#-*- coding: ISO-8859-1 -*-
#pylint: disable-msg=C0301,W0702,C0103,R0903,R0914,R0912,R0915,R0913,C0302,W0141
"""
File       : cmssh_install.py
Author     : Valentin Kuznetsov [ vkuznet AT gmail DOT com ]
Description: cmssh installation script

Some useful URLs
http://www.globus.org/ftppub/gt5/5.0/5.0.4/installers/src/gt5.0.4-all-source-installer.tar.bz2
http://www.nikhef.nl/pub/projects/grid/gridwiki/index.php/Using_voms-proxy-init_on_an_OSX_(10.4_or_higher)_system
https://twiki.grid.iu.edu/bin/view/ReleaseDocumentation/VomsInstallGuide

"""
import sys
if sys.version_info < (2, 6):
    raise Exception("cmssh requires python 2.6 or higher")

# system modules
import os
import re
import stat
import copy
import time
import shutil
import socket
import urllib2
import tarfile
import subprocess

# local modules
from optparse import OptionParser

def osx_test():
    "Test XCode presence for OSX"
    for tool in ['/usr/bin/make', '/usr/include/limits.h']:
        if  not os.path.isfile(tool):
            msg = """######## ERROR ########
Can't locate %s

In order to install and run cmssh on Mac OS X, your system must
have the following components:

Snow Leopard (OSX version 10.8.x)
    - Install Xcode from your installation media (CD/DVD)

Lion (OSX version 10.7.x)
    - Install Xcode 4.2.x from your installation media (CD/DVD) or
      from the Mac App Store

Mountain Lion (OSX version 10.8.x)
    - Install Xcode 4.4.x from the Mac App Store

    - drag the Xcode app to your /Applications folder,
      then run it - It will install almost everything needed

    - Open Xcode preferences to install command line tools, see
      http://bit.ly/TIdOgS

    - Make sure that Java is installed, see http://bit.ly/NoAmkv
""" % tool
            print msg
            sys.exit(1)

def osx_ver():
    "Determine OSX version"
    cmd = 'sw_vers -productVersion'
    res  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    ver = res.stdout.read().replace('\n', '').strip()
    return '.'.join(ver.split('.')[:2])

def get_scram_arch():
    "Determine SCRAM architecture"
    if  os.uname()[0] == 'Darwin':
        osx_test()
        # For OSX Lion I need osx107, but existing packages
        # have problem with bootstrap, see
        # https://hypernews.cern.ch/HyperNews/CMS/get/sw-develtools/1743/1/1/1/1/2/1.html
        # I need to wait until it is fixed, therefore I use arch osx106_amd64_gcc421
        # which seems to be working
        arch = 'osx106_amd64_gcc421'
        if  osx_ver() == '10.7':
            arch = 'osx107_amd64_gcc462'
        if  osx_ver() == '10.8':
            arch = 'osx108_amd64_gcc470'
    elif os.uname()[0] == 'Linux':
        if  os.uname()[-1] == 'x86_64':
            arch = 'slc5_amd64_gcc462'
        else:
            arch = 'slc5_ia32_gcc434'
        # test presence of readline
        if  hasattr(subprocess, "check_output"):
            cond = re.search("libreadline\.so\.5", subprocess.check_output(["/sbin/ldconfig", "-p"]))
        else:
            cmd  = "/sbin/ldconfig -p"
            res  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            libs = [l for l in res.stdout.readlines() if l.find('libreadline.so.5') != -1]
            cond = len(libs)
        if  not cond:
            msg  = 'cmssh on Linux requires readline5. Please verify that'
            msg += ' you have it installed on your system.'
            raise Exception(msg)
    else:
        print 'Unsupported platform'
        sys.exit(1)
    return arch

# generate pkgconfig pc file
def pc_file(path, name, ver, lib, inc):
    "Generate pkgconfig file"
    tmpl = """
prefix=%(path)s
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include/lib%(lib)s

Name: %(name)s
Description: %(name)s
Version: %(ver)s
Libs: -L${libdir} -l%(lib)s
Libs.private: -lz -lbz2
Cflags: -I${includedir} -I${includedir}/%(inc)s
""" % {'path':path, 'name':name, 'lib': lib, 'inc':inc, 'ver': ver}
    return tmpl

# helper function to make natural sort
def try_int(sss):
    "Convert to integer if possible."
    try:
        return int(sss)
    except:
        return sss

def natsort_key(sss):
    "Used internally to get a tuple by which s is sorted."
    return map(try_int, re.findall(r'(\d+|\D+)', sss))

def natcmp(aaa, bbb):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(aaa), natsort_key(bbb))

def natsort(seq, compare=natcmp):
    "In-place natural string sort."
    seq.sort(compare)

def natsorted(seq, compare=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    temp = copy.copy(seq)
    natsort(temp, compare)
    return temp

# Retrieve default scram arch based on OS version and use it globally
# in the rest of the code
DEF_SCRAM_ARCH = get_scram_arch()

# IMPORTANT
# The 2.6.4 version of CMSSW python has bug in OpenSSL
# see https://hypernews.cern.ch/HyperNews/CMS/get/sw-develtools/1667/1.html
# We need to avoid it, otherwise usage of HTTPS will be broken
# So, the osx106_amd64_gcc462 has broken python 2.6.4
# the osx106_amd64_gcc461 has correct python 2.6.4, but broken 2.6.4-cmsX
# the osx106_amd64_gcc421 has corrent python 2.6.4, but it picks root 5.30.02
# which does not have pyROOT library

def find_pkg(pkg, rel_pkgs):
    "Find given package in list of release packages"
    for rel_pkg in rel_pkgs:
        if  pkg == rel_pkg or rel_pkg.find(pkg) != -1:
            return rel_pkg

def find_cms_package(apt_init, pkg, debug=None, lookup=None, rel_pkgs=[]):
    """
    Find latest version of given package in CMSSW repository.
    """
    if  not lookup:
        lookup = pkg
    cmd  = 'source %s; apt-cache search %s | grep "%s" | grep -v toolfile ' % (apt_init, pkg, lookup)
    if  debug:
        print cmd
    if  DEF_SCRAM_ARCH == 'osx107_amd64_gcc462':
        rel_pkg = find_pkg(lookup, rel_pkgs)
        if  rel_pkg:
            name = rel_pkg
        elif  pkg == 'coral':
            # fix coral lib issue for this arch
            # https://hypernews.cern.ch/HyperNews/CMS/get/softwareDistrib/682.html
            name = 'cms+coral+CORAL_2_3_21-cms25'
        else:
            res  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            vers = [r.replace('\n', '').split()[0] for r in res.stdout.readlines()]
            name = natsorted(vers)[-1]
    elif  DEF_SCRAM_ARCH == 'osx106_amd64_gcc421': # Snow Leopard
        if  pkg == 'root':
            name = 'lcg+root+5.30.02-cms4'
        elif pkg == 'py2-matplotlib':
            name = 'external+py2-matplotlib+1.0.1-cms3'
        elif pkg == 'py2-scipy':
            name = 'external+py2-scipy+0.8.0-cms3'
        elif pkg == 'coral':
            name = 'cms+coral+CORAL_2_3_12-cms30'
        elif pkg == 'py2-pycurl':
            name = 'external+py2-pycurl+7.19.0'
        else:
            res  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            vers = [r.replace('\n', '').split()[0] for r in res.stdout.readlines()]
            name = natsorted(vers)[-1]
    else:
        rel_pkg = find_pkg(lookup, rel_pkgs)
        if  rel_pkg:
            name = rel_pkg
        else:
            res  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            vers = [r.replace('\n', '').split()[0] for r in res.stdout.readlines()]
            name = natsorted(vers)[-1]
    return name

def find_installed_pkg(name, debug=0):
    "Find latest version (via natural sort) of installed package for a given name"
    cmd  = 'find $VO_CMS_SW_DIR/$SCRAM_ARCH/%s/*/etc/profile.d -name init.sh' % name
    res  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    script = '/etc/profile.d/init.sh'
    vers = [r.replace('\n', '').split()[0].replace(script, '') for r in res.stdout.readlines()]
    try:
        init = natsorted(vers)[-1] + script
    except:
        if  debug:
            print "Unable to process, name=%s, ver=%s" % (name, vers)
        if  vers:
            print natsorted(vers)
        return None, None, None
    root = '/'.join(init.split('/')[:-3])
    ver  = init.split('/')[-4]
    return init, root, ver

def available_architectures():
    "Fetch CMSSW drivers"
    arch = os.uname()[0]
    pat1 = re.compile('.*-driver.txt</a>.*')
    pat2 = re.compile('^[osx,slc].*')
    url  = 'http://cmsrep.cern.ch/cmssw/cms/'
    data = urllib2.urlopen(url)
    for line in data.readlines():
        if  pat1.match(line):
            line = line.split('</a>')[0].split('">')[-1]
            if  pat2.match(line):
                if  arch == 'Linux' and line[:3] == 'slc':
                    yield line.replace('-driver.txt', '')
                elif arch == 'Darwin' and line[:3] == 'osx':
                    yield line.replace('-driver.txt', '')

def install_pip_pkg(pkgs, cms_env, path, debug, pkg, ver=None, opts=None, env_list=None):
    """
    Install pip package. User must provide:
        - pkgs, installed package list
        - cms_env, cms_environment source string
        - path, cmssh install path
        - debug flag
        - pkg to be installed
        - optional: version, options and environment list (key-value tuple)
    """
    print "Install", pkg
    cmd = cms_env
    if  env_list:
        for key, val in env_list:
            os.environ[key] = val
            if  debug:
                print "os.environ[%s]=%s" % (key, val)
    if  not pkgs.has_key(pkg):
        if  ver:
            cmd += '%s/install/bin/pip install %s==%s' % (path, pkg, ver)
        else:
            cmd += '%s/install/bin/pip install --upgrade %s' % (path, pkg)
        if  opts:
            cmd += ' ' + opts
        exe_cmd(path, cmd, debug, log='%s.log' % pkg)
    if  env_list:
        for key, val in env_list:
            del os.environ[key] # clean-up env settings

def replace_in_file(fname, pat, new_pat):
    "Replace given pattern to new one for given file name"
    new_fname = '%s.%s' % (fname, int(time.time()))
    with open(fname, 'r') as old_file:
        with open(new_fname, 'w') as new_file:
            lines = old_file.read()
            lines = lines.replace(pat, new_pat)
            new_file.write(lines)
    os.remove(fname)
    shutil.move(new_fname, fname)

def siteconfig():
    "Generate dummy site config file used by lumiDB"
    header = """<site-local-config><site name="T3_XX_YYYY"><calib-data><frontier-connect>"""
    bottom = """</frontier-connect><catalog url=""/></calib-data></site></site-local-config>"""
    body   = ""
    idict  = {'proxy': {'host': 'cmst0frontier', 'port': 3128},
              'backupproxy': {'host': 'cmsbpfrontier', 'port': 3128},
              'server': {'host': 'cmsfrontier', 'port': 8000, 'api': 'FrontierInt'}}
    for tag, item in idict.items():
        for num in ['', '1', '2']:
            url = "http://%s%s.cern.ch:%s" % (item['host'], num, item['port'])
            if  item.has_key('api'):
                url += "/%s" % item['api']
            entry = '<%s url="%s"/>' % (tag, url)
            body += entry
    return header + body + bottom

class MyOptionParser:
    """option parser"""
    def __init__(self):
        self.parser = OptionParser()
        self.parser.add_option("-v", "--verbose", action="store",
            type="int", default=0, dest="debug",
            help="verbose output")
        self.parser.add_option("-d", "--dir", action="store",
            type="string", default=None,
            dest="install_dir", help="install directory")
        self.parser.add_option("-i", "--install", action="store_true",
            dest="install", help="install cmssh and its dependencies")
        self.parser.add_option("-u", "--upgrade", action="store_true",
            dest="upgrade", help="upgrade cmssh")
        self.parser.add_option("--version", action="store",
            type="string", default="1.0.0",
            dest="version", help="get specific version of cmssh, e.g. master")
        drivers = ', '.join(available_architectures())
        self.parser.add_option("--arch", action="store",
            type="string", default=None, dest="arch",
            help="CMSSW architectures:\n%s, default %s" % (drivers, DEF_SCRAM_ARCH))
        self.parser.add_option("--cmssw", action="store",
            type="string", default=None, dest="cmssw",
            help="specify location of CMSSW install area")
        if  socket.gethostname().find('cern.ch') != -1:
            lcg_default = '/afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh'
        else:
            lcg_default = ''
        self.parser.add_option("--lcg", action="store",
            type="string", default=lcg_default, dest="lcg",
            help="specify script to setup LCG/OSG environment")
        self.parser.add_option("--multi-user", action="store_true",
            default=False, dest="multi_user",
            help="install cmssh in multi-user environment")
        self.parser.add_option("--unsupported", action="store_true",
            dest="unsupported",
            help="enforce installation on unsupported CMS platforms, e.g. non SLC Linux")
        self.parser.add_option("--seed", action="store",
            type="string", default="", dest="seed",
            help="seed dependencies from given release")

    def get_opt(self):
        """Returns parse list of options"""
        return self.parser.parse_args()

def getdata(url, verbose=0):
    """Invoke URL call and retrieve data for given url/params/headers"""
    headers = {}
    if  verbose:
        print '+++ getdata url=%s' % url
    req = urllib2.Request(url)
    for key, val in headers.items():
        req.add_header(key, val)
    if  verbose > 1:
        handler = urllib2.HTTPHandler(debuglevel=1)
        opener  = urllib2.build_opener(handler)
        urllib2.install_opener(opener)
    data = urllib2.urlopen(req)
    return data.read()

def is_installed(url, path):
    "Check if already installed a package for given URL"
    fname = os.path.join(path, '.packages')
    if  os.path.isfile(fname):
        with open(fname, 'r') as packages:
            for line in packages.readlines():
                if  url == line.replace('\n', ''):
                    return True
    return False

def add_url2packages(url, path):
    "Add url to packages file"
    with open(os.path.join(path, '.packages'), 'a') as packages:
        packages.write(url + '\n')

def get_file(url, fname, path, debug, ext='r:gz'):
    """Fetch tarball from given url and store it as fname, untar it into given path"""
    os.chdir(path)
    with open(fname, 'w') as tar_file:
        tar_file.write(getdata(url, debug))
    tar = tarfile.open(fname, ext)
    top_names = set([r.split('/')[0] for r in tar.getnames()])
    if  len(top_names) == 1:
        dir_name = top_names.pop()
        if  os.path.isdir(dir_name):
            try:
                os.removedirs(dir_name)
            except:
                pass
    tar.extractall(path)
    tar.close()
    add_url2packages(url, path)

def exe_cmd(idir, cmd, debug, msg=None, log='install.log'):
    """Execute given command in a given dir"""
    if  msg:
        print msg
    os.chdir(idir)
    if  debug:
        print "cd %s\n%s" % (idir, cmd)
    with open(log, 'w') as logstream:
        try:
            retcode = subprocess.call(cmd, shell=True, stdout=logstream, stderr=logstream)
            if  retcode < 0:
                print >> sys.stderr, "Child was terminated by signal", -retcode
        except OSError, err:
            print >> sys.stderr, "Execution failed:", err

def check_system(unsupported):
    "Check system requirements"
    # check presence of Java, required for GRID middleware
    if  not os.environ.has_key('JAVA_HOME'):
        print "JAVA_HOME environment is required to install GRID middleware tools"
        print "Please install Java and appropriately setup JAVA_HOME"
        print "For example, export JAVA_HOME=/usr"
        sys.exit(1)

    # check Ubuntu default shell
    if  os.uname()[3].find('Ubuntu') != -1 or unsupported:
        if  os.readlink('/bin/sh') != 'bash':
            msg  = 'The /bin/sh is pointing to %s.\n'
            msg += 'For proper installation of CMSSW software\n'
            msg += 're-configure /bin/sh to point to /bin/bash.\n'
            msg += 'On Ubuntu, if you have /bin/sh -> /bin/dash, just do:\n'
            msg += 'sudo dpkg-reconfigure dash'
            print msg
            sys.exit(1)

def test_Fortran():
    "Test presence of Fortran on the system"
    for cmd in ['g95', 'f95', 'f90', 'f77', \
            'xlf90', 'xlf', 'ifort', 'ifc', \
            'g77', 'gfortran', 'pgfortran']:
        cmd = 'which %s' % cmd
        res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        rpath = res.stdout.read()
        return rpath

def main():
    "Main function"
    mgr = MyOptionParser()
    opts, _ = mgr.get_opt()
    use_lcg = opts.lcg

    platform = os.uname()[0]
    if  platform == 'Darwin':
        if  not os.environ.has_key('JAVA_HOME'):
            os.environ['JAVA_HOME'] = '/Library/Java/Home'
    elif platform == 'Linux':
        if  not os.environ.has_key('JAVA_HOME'):
            if  os.path.isfile('/usr/bin/java') or os.path.islink('/usr/bin/java'):
                os.environ['JAVA_HOME'] = '/usr'
            elif os.path.isfile('/usr/local/bin/java') or os.path.islink('/usr/local/bin/java'):
                os.environ['JAVA_HOME'] = '/usr/local'
            elif os.path.isfile('/opt/local/bin/java') or os.path.islink('/opt/local/bin/java'):
                os.environ['JAVA_HOME'] = '/opt/local'

    if  not opts.install and not opts.upgrade:
        print "Usage: cmssh_install.py --help"
        sys.exit(0)
    check_system(opts.unsupported)
    debug = opts.debug
    idir = opts.install_dir
    if  not idir:
        msg  = "Please specify the install area"
        msg += " (it should have enough space to hold CMSSW releases)"
        print msg
        sys.exit(1)
    if  idir[0] != '/':
        msg = "Please specify absolute path for install area, e.g. /path/%s" % idir
        print msg
        sys.exit(1)
    arch   = opts.arch
    path   = os.path.join(idir, 'soft')

    # setup platform and unsupported flag
    unsupported_linux = False
    if  os.uname()[3].find('Ubuntu') != -1 or opts.unsupported:
        unsupported_linux = True

    # setup system architecture
    parch = 'x86'
    arch  = opts.__dict__.get('arch', None)
    if  platform == 'Linux':
        parch = 'x86_64'
        if  unsupported_linux:
            vdt_ver = 'deb_5.0'
        else:
            vdt_ver = 'rhap_5'
        if  not arch:
            arch = DEF_SCRAM_ARCH
    elif platform == 'Darwin':
        vdt_ver  = 'macos_10.4'
        if  not arch:
            arch = DEF_SCRAM_ARCH
    else:
        print 'Unsupported OS "%s"' % platform
        sys.exit(1)
    if  not arch:
        print "Unsupported architecture"
        sys.exit(1)

    print 'Checking CMSSW'
    if  debug:
        print 'Probe architecture', arch
    try:
        os.makedirs(os.path.join(path, 'logs'))
    except:
        pass
    os.chdir(path)
    os.environ['LANG'] = 'C'
    sdir = '%s/CMSSW' % path

    if  opts.cmssw:
        # check if default architecture is present
        if  arch in os.listdir(opts.cmssw):
            if  not os.path.islink(sdir):
                os.symlink(opts.cmssw, sdir)
            os.environ['SCRAM_ARCH'] = arch
            os.environ['VO_CMS_SW_DIR'] = sdir
            if  debug:
                print 'Will use %s/%s' % (sdir, arch)
        else:
            print "Please supply via --arch the architecture from %s you wish to use"\
                % opts.cmssh
            sys.exit(1)
    else: # do local CMSSW bootstrap
        try:
            os.makedirs(sdir)
        except:
            pass
        os.environ['VO_CMS_SW_DIR'] = sdir
        os.environ['SCRAM_ARCH'] = arch
        url  = 'http://cmsrep.cern.ch/cmssw/cms/bootstrap.sh'
        os.chdir(path)
        if  not is_installed(url, path):
            os.chdir(sdir)
            with open('bootstrap.sh', 'w') as bootstrap:
                bootstrap.write(getdata(url, debug))
            if  os.uname()[0].lower() == 'linux':
                os.rename('bootstrap.sh', 'b.sh')
                cmd = 'cat b.sh | sed "s,\$seed \$unsupportedSeeds,\$seed \$unsupportedSeeds libreadline5,g" > bootstrap.sh'
                subprocess.call(cmd, shell=True)
            os.chmod('bootstrap.sh', 0755)
            cmd  = 'sh -x $VO_CMS_SW_DIR/bootstrap.sh setup -path $VO_CMS_SW_DIR -arch $SCRAM_ARCH'
            if  unsupported_linux:
                cmd += ' -unsupported_distribution_hack'
            os.chdir(path)
            exe_cmd(sdir, cmd, debug, 'Bootstrap CMSSW', log='bootstrap.log')
            apt_init, _root, _ver = find_installed_pkg('external/apt', debug)
            cmd  = 'source %s' % apt_init
            cmd += 'apt-get install external+fakesystem+1.0; '
            cmd += 'apt-get update; '
            exe_cmd(sdir, cmd, debug, 'Init CMSSW apt repository', log='aptget.log')
            # get latest non-patched/non-pre release and find out package versions
            latest_release = opts.seed
            if  not latest_release:
                cmd  = "source %s; apt-cache search CMSSW_ " % apt_init
                cmd += '| egrep -v -i "fwlite|pre|patch|dqm"'
                cmd += "| awk '{print $1}'"
                cmd += '| grep "^cms[+]cmssw[+]CMSSW_[0-9]_[0-9]_[0-9]$"'
                cmd += "| tail -1"
                exe_cmd(sdir, cmd, debug, 'Find latest CMSSW release', log='cmssw_rel.log')
                with open('cmssw_rel.log', 'r') as stream:
                    latest_release = stream.read().replace('\n', '').strip()
            rel_pkgs = []
            if  latest_release:
                cmd  = 'source %s; echo "N" | apt-get install %s' % (apt_init, latest_release)
                exe_cmd(sdir, cmd, debug, 'Find %s release deps' % latest_release, log='cmssw_rel_dep.log')
                with open('cmssw_rel_dep.log', 'r') as stream:
                    rel_pkgs = [p for p in stream.read().split() if p.find('+')!=-1]
                if  debug:
                    print "Packages:"
                    for pkg in rel_pkgs:
                        print pkg
            # install useful set of CMS libraries
            cms_libs = ['root', 'coral', 'py2-pycurl', 'py2-matplotlib', 'py2-scipy']
            if  platform == 'Darwin' and osx_ver() == '10.6':
                cms_libs += ['freetype']
            for cmspkg in cms_libs:
                msg  = 'Install CMSSW %s' % cmspkg
                if  cmspkg == 'root':
                    lookup = 'lcg+root'
                elif  cmspkg == 'coral':
                    lookup = 'coral+CORAL'
                else:
                    lookup = ''
                name = find_cms_package(apt_init, cmspkg, debug, lookup, rel_pkgs)
                cmd  = 'source %s; echo "Y" | apt-get install %s' % (apt_init, name)
                log  = '%s/logs/%s.log' % (path, cmspkg)
                exe_cmd(sdir, cmd, debug, msg, log)
            # add bootstrap url into soft/.packages
            add_url2packages(url, path)

    # command to setup CMSSW python
    _coral_init, coral_root, coral_ver = find_installed_pkg('cms/coral', debug)
    python_init, python_root, python_ver = find_installed_pkg('external/python', debug)
    if  not python_init.find('init.sh') != -1:
        msg  = '\nUnable to locate python in:'
        msg += '\n%s/%s/external/python' % (os.environ['VO_CMS_SW_DIR'], os.environ['SCRAM_ARCH'])
        msg += '\nPlease check CMSSW area and/or choose another architecture\n'
        print msg
        sys.exit(1)
    pver     = '.'.join(python_ver.split('.')[0:2])
    install_dir = '%s/install/lib/python%s/site-packages' % (path, pver)
    os.environ['PYTHONPATH'] = install_dir
    try:
        os.makedirs(install_dir)
    except:
        pass
    cms_env  = 'source %s;' % python_init
    if  debug:
        print "CMSSW python: %s/%s" % (python_root, python_ver)
        print "python version", pver

    if  platform == 'Darwin' and osx_ver() == '10.6':
        # CMSSW pcre is too old and srm software uses grep which linked to newer
        # pcre library, therefore install pcre 7.9 which is suitable for this case
        print "Install pcre"
        ver = '7.9'
        url = 'http://downloads.sourceforge.net/pcre/%s/pcre-%s.tar.gz' % (ver, ver)
        if  not is_installed(url, path):
            get_file(url, 'expat-%s.tar.gz' % ver, path, debug)
            cmd = cms_env + './configure --prefix=%s/install; make; make install' % path
            os.chdir(os.path.join(path, 'pcre-%s' % ver))
            exe_cmd(os.path.join(path, 'pcre-%s' % ver), cmd, debug, log='pcre.log')

    print "Install expat"
    ver = '2.0.1'
    url = 'http://sourceforge.net/projects/expat/files/expat/%s/expat-%s.tar.gz/download?use_mirror=iweb' % (ver, ver)
    if  not is_installed(url, path):
        get_file(url, 'expat-%s.tar.gz' % ver, path, debug)
        if  parch == 'x86':
            cflags = 'CFLAGS=-m32'
        else:
            cflags = ''
        cmd = cms_env + '%s ./configure --prefix=%s/install; make; make install' % (cflags, path)
        os.chdir(os.path.join(path, 'expat-%s' % ver))
        exe_cmd(os.path.join(path, 'expat-%s' % ver), cmd, debug, log='expat.log')

    print "Install PythonUtilities"
    url = "http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/FWCore/PythonUtilities.tar.gz?view=tar"
    if  not is_installed(url, path):
        get_file(url, 'PythonUtilities.tar.gz', path, debug)
        cmd = 'touch __init__.py; mv python/*.py .'
        exe_cmd(os.path.join(path, 'PythonUtilities'), cmd, debug)
        os.chdir(path)
        cmd = 'mkdir FWCore; touch FWCore/__init__.py; mv PythonUtilities FWCore'
        exe_cmd(path, cmd, debug)

#    print "Install CRAB3"
#    ver = '3.0.6a'
#    url = 'http://cmsrep.cern.ch/cmssw/comp/SOURCES/slc5_amd64_gcc461/cms/crab-client3/%s/crabclient3.tar.gz' % ver
#    if  not is_installed(url, path):
#        get_file(url, 'crabclient3.tar.gz', path, debug)

    print "Install CRAB"
    os.chdir(path)
    crab_ver = 'CRAB_current'
    url = 'http://cmsdoc.cern.ch/cms/ccs/wm/www/Crab/Docs/%s.tgz' % crab_ver
    if  not is_installed(url, path):
        try:
            get_file(url, 'crab.tar.gz', path, debug)
            cmd = 'cd %s; ./configure' % crab_ver
            exe_cmd(path, cmd, debug, log='crab.log')
        except:
            print "Unable to fetch CRAB tar ball, will skip this step"

    print "Install WMCore"
    ver = '0.8.21'
    url = 'http://cmsrep.cern.ch/cmssw/comp/SOURCES/slc5_amd64_gcc461/cms/wmcore/%s/WMCORE.tar.gz' % ver
    if  not is_installed(url, path):
        get_file(url, 'wmcore.tar.gz', path, debug)

    print "Install certificates"
    url = 'http://vdt.cs.wisc.edu/software/certificates/62/certificates-62-1.tar.gz'
    url = 'http://dist.eugridpma.info/distribution/igtf/current/accredited/igtf-preinstalled-bundle-classic.tar.gz'
    opath = os.getcwd()
    cpath = os.path.join(path, 'certificates')
    if  not is_installed(url, cpath):
        os.makedirs(cpath)
        get_file(url, 'certificates.tar.gz', cpath, debug)
    os.chdir(opath)

    # test local setup of GRID middleware (LCG/OSG)
    if  platform == 'Linux' and use_lcg:
        print "Skip GRID middleware install, use local setup"
    else:
        print "Install Globus"
        os.chdir(path)
        url = 'http://vdt.cs.wisc.edu/software/globus/4.0.8_VDT2.0.0gt4nbs/vdt_globus_essentials-VDT2.0.0-3-%s_%s.tar.gz' % (parch, vdt_ver)
        if  not is_installed(url, path):
            get_file(url, 'globus.tar.gz', path, debug)

        print "Install Myproxy"
        url = 'http://vdt.cs.wisc.edu/software/myproxy/5.3_VDT-2.0.0/myproxy_client-5.3-%s_%s.tar.gz' % (parch, vdt_ver)
        if  not is_installed(url, path):
            get_file(url, 'myproxy_client.tar.gz', path, debug)
        url = 'http://vdt.cs.wisc.edu/software/myproxy/5.3_VDT-2.0.0/myproxy_essentials-5.3-%s_%s.tar.gz' % (parch, vdt_ver)
        if  not is_installed(url, path):
            get_file(url, 'myproxy_essentials.tar.gz', path, debug)

        print "Install VOMS"
        url = 'http://vdt.cs.wisc.edu/software/voms/1.8.8-2p1-1/voms-client-1.8.8-2p1-%s_%s.tar.gz' % (parch, vdt_ver)
        if  not is_installed(url, path):
            get_file(url, 'voms-client.tar.gz', path, debug)
        url = 'http://vdt.cs.wisc.edu/software/voms/1.8.8-2p1-1/voms-essentials-1.8.8-2p1-%s_%s.tar.gz' % (parch, vdt_ver)
        if  not is_installed(url, path):
            get_file(url, 'voms-essentials.tar.gz', path, debug)

        print "Install LCG info"
        url = 'http://vdt.cs.wisc.edu/software/lcg-infosites/2.6-2/lcg-infosites-2.6-2.tar.gz'
        if  not is_installed(url, path):
            get_file(url, 'lcg-infosites.tar.gz', path, debug)
        url = 'http://vdt.cs.wisc.edu/software/lcg-info//1.11.4-1/lcg-info-1.11.4-1.tar.gz'
        if  not is_installed(url, path):
            get_file(url, 'lcg-info.tar.gz', path, debug)

        print "Install SRM client"
        ver = '2.2.1.3.19'
        url = 'http://vdt.cs.wisc.edu/software/srm-client-lbnl/%s/srmclient2-%s.tar.gz' \
            % (ver, ver)
        if  not is_installed(url, path):
            get_file(url, 'srmclient.tar.gz', path, debug)
            cmd  = cms_env + './configure --with-java-home=$JAVA_HOME --enable-clientonly'
            cmd += ' --with-globus-location=%s/globus' % path
            cmd += ' --with-cacert-path=%s/certificates' % path
            cmd += ' --with-srm-home=%s/srmclient2' % path
            exe_cmd(os.path.join(path, 'srmclient2/setup'), cmd, debug, log='srmclient.log')
            # fix Lion, Unable to load native library: libjava.jnilib problem
            # http://stackoverflow.com/questions/1482450/broken-java-mac-10-6
            if  platform == 'Darwin' and osx_ver() == '10.6':
                pat = 'export CLASSPATH'
                new_pat = pat + '\nunset DYLD_LIBRARY_PATH\n'
                for fname in os.listdir(os.path.join(path, 'srmclient2/bin')):
                    if  fname.find('srm-') != -1:
                        filename = os.path.join(path, 'srmclient2/bin/' + fname)
                        replace_in_file(filename, pat, new_pat)
                        os.chmod(filename, 0755)

    print "Install zmq"
    os.chdir(path)
    zmq_ver = '2.2.0'
    url = 'http://download.zeromq.org/zeromq-%s.tar.gz' % zmq_ver
    if  not is_installed(url, path):
        get_file(url, 'zmq.tar.gz', path, debug) # it call add_url2packages
        cmd = 'cd zeromq-%s; ./configure --prefix=%s/install' % (zmq_ver, path)
        cmd += '; make install'
        exe_cmd(path, cmd, debug, log='zmq.log')

    print "Install setuptools"
    os.chdir(path)
    s_ver = '0.6c11'
    url = 'http://pypi.python.org/packages/source/s/setuptools/setuptools-%s.tar.gz' % s_ver
    if  not is_installed(url, path):
        get_file(url, 'setuptools.tar.gz', path, debug) # it call add_url2packages
        cmd = cms_env + 'cd setuptools-%s; python setup.py install --prefix=%s/install' % (s_ver, path)
        exe_cmd(path, cmd, debug, log='setuptools.log')

    print "Install pip"
    os.chdir(path)
    ver = '1.10.1'
    url = 'https://pypi.python.org/packages/source/v/virtualenv/virtualenv-%s.tar.gz' % ver
    if  not is_installed(url, path):
        get_file(url, 'virtualenv.tar.gz', path, debug) # it call add_url2packages
        cmd = cms_env + 'python virtualenv.py %s/install' % path
        cmd += '; . %s/install/activate' % path
        idir = '%s/virtualenv-%s' % (path, ver)
        exe_cmd(idir, cmd, debug, log='pip.log')
        add_url2packages(url, path)
#    url = 'https://raw.github.com/pypa/virtualenv/master/virtualenv.py'
#    if  not is_installed(url, path):
#        with open('virtualenv.py', 'w') as fname:
#            fname.write(getdata(url, debug))
#        cmd = cms_env + 'python %s/virtualenv.py %s/install' % (path, path)
#        cmd += '; . %s/install/activate' % path
#        exe_cmd(path, cmd, debug, log='pip.log')
#        add_url2packages(url, path)

    # get list of installed packages in pip repository
    cmd = cms_env + '%s/install/bin/pip freeze' % path
    exe_cmd(path, cmd, debug, log='installed.pkg')
    pip_packages = {}
    logfile = '%s/installed.pkg' % path
    with open(logfile, 'r') as installed_packages:
        for line in installed_packages.readlines():
            if  line.lower().find('warning') != -1 or line.find('#') != -1:
                continue
            if  not line:
                continue
            try:
                pkg, ver = line.replace('\n', '').split('==')
                pip_packages[pkg] = ver
            except Exception as exp:
                print "Fail at line '%s'" % line
                print exp
                raise

    pkg = 'ipython'
    install_pip_pkg(pip_packages, cms_env, path, debug, pkg)
    # fix pylab message
    fname = '%s/install/lib/python%s/site-packages/IPython/core/pylabtools.py' \
        % (path, pver)
    content = None
    with open(fname, 'r') as source:
        content = source.read()
        content = content.replace('Welcome to pylab, a matplotlib-based', 'cmssh+pylab')
        content = content.replace("For more information, type 'help(pylab)'.", '')
    if  content:
        with open(fname, 'w') as output:
            output.write(content)

    # update local pkgconfig with libpng/freetype
    try:
        os.makedirs('%s/install/lib/pkgconfig' % path)
    except:
        pass
    with open('%s/install/lib/pkgconfig/libpng.pc' % path, 'w') as stream:
        png_init, png_root, png_ver = find_installed_pkg('external/libpng', debug)
        prefix = '%s/install' % path
        name   = 'libpng'
        lib    = 'png%s' % ''.join(png_ver.split('.')[:2])
        ver    = png_ver
        inc    = ''
        stream.write(pc_file(prefix, name, ver, lib, inc))
    with open('%s/install/lib/pkgconfig/freetype2.pc' % path, 'w') as stream:
        ft_init, ft_root, ft_ver = find_installed_pkg('external/freetype', debug)
        prefix = '%s/install' % path
        name   = 'freetype'
        lib    = 'freetype'
        ver    = ft_ver
        inc    = 'freetype2'
        stream.write(pc_file(prefix, name, ver, lib, inc))
    # install standard libraries
    std_pkgs = ['Routes', 'python-dateutil', 'decorator',
            'pyOpenSSL', 'paramiko', 'pyzmq', 'tornado', 'pytz', 'pycrypto',
            'numpy', 'matplotlib', 'html2text', 'feedparser', 'jinja2',
    ]
    for pkg in std_pkgs:
        ver  = None
        args = None
        env_list = []
        if  pkg.lower() == 'feedparser':
            ver = '5.1.2' # tested
        if  pkg.lower() == 'pyopenssl':
            # use 0.12 version of pyOpenSSL due to
            # http://stackoverflow.com/questions/7340784/easy-install-pyopenssl-error
            ver = '0.12'
        if  pkg.lower() == 'pyzqm':
            # add intstall path option for zmq
            args = '--install-option="--zmq=%s/install"' % path
        if pkg == 'matplotlib':
            try:
                if  platform == 'Darwin' and \
                    os.path.isfile('/usr/bin/llvm-gcc') and \
                    os.path.isfile('/usr/bin/llvm-g++'):
                    env_list  = [('CC', 'llvm-gcc'), ('CXX', 'llvm-g++')]
                    env_list += [('MACOSX_DEPLOYMENT_TARGET', osx_ver())]
                    ft_init, ft_root, _   = find_installed_pkg('external/freetype', debug)
                    png_init, png_root, _ = find_installed_pkg('external/libpng', debug)
                    cms_env2 = 'source %s; source %s; source %s;' \
                            % (python_init, ft_init, png_init)
                    cflags  = '-I%s/include -I%s/include/freetype2 -I%s/include' \
                            % (ft_root, ft_root, png_root)
                    ldflags = '-L%s/lib -L%s/lib' % (ft_root, png_root)
                    env_list += [('CFLAGS', cflags), ('LDFLAGS', ldflags)]
                    pkgconfig = '%s/install/lib/pkgconfig' % path
                    env_list += [('PKG_CONFIG_PATH', pkgconfig)]
                    install_pip_pkg(pip_packages, cms_env2, path, debug, pkg, ver, args, env_list)
            except:
                pass
        else:
            install_pip_pkg(pip_packages, cms_env, path, debug, pkg, ver, args, env_list)

    # install readline after pip, since it requires setuptools
    print "Install readline"
#    if  platform == 'Darwin' and not is_installed(url, path):
#        if  pver == '2.7':
#            url = 'http://pypi.python.org/packages/2.7/r/readline/readline-6.2.2-py2.7-macosx-10.7-intel.egg'
#            md5 = '25383d860632d4a1521961ba68a52fe2'
#        if  pver == '2.6':
#            url = 'http://pypi.python.org/packages/2.6/r/readline/readline-6.2.2-py2.6-macosx-10.6-universal.egg'
#            md5 = '13d63c76be4ff09c5d55fbbe3b6ab2c7'
#        with open('readline.egg', 'w') as readline:
#            readline.write(getdata(url, debug))
#        cmd = 'cd %s/install/lib/python%s/site-packages; unzip -n %s/readline.egg' % (path, pver, path)
#        exe_cmd(path, cmd, debug, log='readline.log')

    ver = '6.2.2'
    url = 'http://pypi.python.org/packages/source/r/readline/readline-%s.tar.gz' % ver
    if  platform == 'Darwin' and not is_installed(url, path):
        get_file(url, 'readline.tar.gz', path, debug)
        cmd = """#!/bin/bash
export CMSSH_ROOT={path}
export VO_CMS_SW_DIR=$CMSSH_ROOT/CMSSW
export SCRAM_ARCH={arch}
export LANG="C"
source {python_init}
idir={path}/install
mkdir -p $idir/lib/python{pver}/site-packages
export PYTHONPATH=$PYTHONPATH:$CMSSH_ROOT/install/lib/python{pver}/site-packages:$idir/lib/python{pver}/site-packages
python setup.py install --prefix=$idir
export CFLAGS='-arch x86_64'
export LDFLAGS='-arch x86_64'
cd readline
./configure CPPFLAGS='-DNEED_EXTERN_PC -fPIC'
make
cd -
python setup.py install --prefix=$idir
""".format(path=path, arch=arch, pver=pver, python_init=python_init)
        exe_cmd(os.path.join(path, 'readline-%s' % ver), cmd, debug, log='readline.log')

    print "Install LumiDB"
    os.chdir(path)
    url = 'http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoLuminosity.tar.gz?view=tar'
    if  not is_installed(url, path):
        get_file(url, 'lumidb.tar.gz', path, debug)
        dst = os.path.join(path, 'install/lib/python%s/site-packages/RecoLuminosity' % pver)
        try:
            os.makedirs(dst)
        except:
            pass
        shutil.copytree('RecoLuminosity/LumiDB/python', os.path.join(dst, 'LumiDB'))
        with open(os.path.join(dst, '__init__.py'), 'w') as init_file:
            init_file.write("")
        shutil.copy(os.path.join(dst, '__init__.py'), os.path.join(dst, 'LumiDB'))

    if  opts.upgrade:
        print "Upgrade cmssh"
    else:
        print "Install cmssh"
    os.chdir(path)
    url = 'http://github.com/dmwm/cmssh/tarball/%s/' % opts.version
    cmssh_ver = [i for i in url.split('/') if i][-1]
    cmssh_ts  = time.strftime("%Y-%m-%d %H:%M:%S GMT", time.gmtime())
    if  opts.upgrade or not is_installed(url, path):
        try:
            cmd = 'rm -rf dmwm-cmssh*; rm -rf cmssh'
            exe_cmd(path, cmd, debug)
        except:
            pass
        try:
            cmd = 'rm -rf .ipython'
            exe_cmd(path, cmd, debug)
        except:
            pass
        get_file(url, 'cmssh.tar.gz', path, debug)
        cmd = 'mv dmwm-cmssh* %s/cmssh' % path
        exe_cmd(path, cmd, debug)

    print "Create configuration"
    os.chdir(path)
    with open('setup.sh', 'w') as setup:
        msg  = '#!/bin/bash\nexport CMSSH_ROOT=%s\n' % path
        msg += """
# Clean-up release environment
dir=$CMSSH_ROOT/install/lib/release_lib
if [ -d $dir ] || [ -L $dir ]; then
   rm -rf $dir
fi
dir=$CMSSH_ROOT/install/lib/release_root
if [ -d $dir ] || [ -L $dir ]; then
   rm -rf $dir
fi
dir=$CMSSH_ROOT/install/lib/release_external
if [ -d $dir ] || [ -L $dir ]; then
   rm -rf $dir
fi
"""
        msg += 'export CMSSH_INSTALL_DIR=$CMSSH_ROOT/install/lib/python%s/site-packages\n' % pver
        msg += 'echo -n "Loading dependencies:"\n'
        msg += """cms_init()
{
    if [ -f $VO_CMS_SW_DIR/$SCRAM_ARCH/$1/$2/etc/profile.d/init.sh ]; then
        echo -n "."
        source $VO_CMS_SW_DIR/$SCRAM_ARCH/$1/$2/etc/profile.d/init.sh
    fi
}
set_cmd()
{
    if  [ "$1" == "lcg" ]; then
        export LCG_CP=`command -v lcg-cp`
        export LCG_LS=`command -v lcg-ls`
        export LCG_RM=`command -v lcg-del`
    fi
    if [ "$1" == "srm" ]; then
        if  [ -n "`command -v srm-ls`" ]; then
            export SRM_CP=`command -v srm-copy`
            export SRM_LS=`command -v srm-ls`
            export SRM_RM=`command -v srm-rm`
            export SRM_MKDIR=`command -v srm-mkdir`
            export SRM_RMDIR=`command -v srm-rmdir`
        fi
        if [ -n "`command -v srmls`" ]; then
            export SRM_CP=`command -v srmcp`
            export SRM_LS=`command -v srmls`
            export SRM_RM=`command -v srmrm`
            export SRM_MKDIR=`command -v srmmkdir`
            export SRM_RMDIR=`command -v srmrmdir`
        fi
    fi
}
link_root()
{
    dir=$CMSSH_ROOT/install/lib/release_root
    if [ -d $dir ] || [ -L $dir ]; then
        rm -rf $dir
    fi
    ln -s $VO_CMS_SW_DIR/$SCRAM_ARCH/$1/$2 $dir
}
coral_init()
{
    coral=$CORAL_DIR/$SCRAM_ARCH
    coral_external=$CORAL_DIR/external/$SCRAM_ARCH/lib
    export PYTHONPATH=$PYTHONPATH:$CMSSH_ROOT:$coral/python:$coral/lib:$coral_external
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$coral/lib:$coral_external
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$coral/lib:$coral_external
    if [ -f $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_PyCoral.dylib ]; then
        if [ ! -f $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_PyCoral.so ]; then
            ln -s $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_PyCoral.dylib $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_PyCoral.so
        fi
    fi
    if [ -f $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_ConnectionService.dylib ]; then
        if [ ! -f $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_ConnectionService.so ]; then
            ln -s $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_ConnectionService.dylib $CORAL_DIR/$SCRAM_ARCH/lib/liblcg_ConnectionService.so
        fi
    fi
}
\n
"""
        msg += 'export VO_CMS_SW_DIR=$CMSSH_ROOT/CMSSW\n'
        msg += 'export SCRAM_ARCH=%s\n' % arch
        msg += 'export LANG="C"\n'
        if  not opts.multi_user:
            msg += 'export CMSSW_RELEASES=$CMSSH_ROOT/Releases\n'
        msg += 'if [ -f $VO_CMS_SW_DIR/cmsset_default.sh ]; then\n'
        msg += '   source $VO_CMS_SW_DIR/cmsset_default.sh\nfi\n'
        msg += "network=`hostname -d 2> /dev/null`\n"
        msg += 'if [ X$network == X"cern.ch" ]; then\n'
        msg += "source /afs/cern.ch/project/eos/installation/pro/etc/setup.sh\n"
        msg += "export CMSSH_EOS=1\nfi\n"
        msg += 'export CRAB_ROOT=$CMSSH_ROOT/%s\n' % crab_ver
        if  not opts.multi_user:
            msg += 'export OLD_PATH=$PATH\n'
            msg += 'export PATH=/usr/bin:/bin:/usr/sbin:/sbin\n'
        msg += 'unset PYTHONPATH\n'
        if  use_lcg:
            msg += '# use local LCG installation\n'
            msg += '. %s\n\n' % use_lcg
        if  os.environ.has_key('LD_LIBRARY_PATH'):
            msg += 'export LD_LIBRARY_PATH=%s\n' % os.environ['LD_LIBRARY_PATH']
        deps = ['external/apt', 'lcg/root', 'external/python',
                'external/xz', 'external/pcre',
                'external/freetype', 'external/libpng', 'external/lapack',
                'external/libjpg', 'external/libtiff', 'external/libungif',
                'external/py2-matplotlib', 'external/py2-scipy',
                'external/py2-numpy', 'external/curl', 'external/py2-pycurl',
                'external/xrootd', 'external/boost', 'cms/coral']
        if  platform == 'Darwin' and osx_ver() == '10.6':
            deps += ['external/xerces-c', 'external/frontier_client']
        matplotlib_ver = None
        for pkg in deps:
            func = 'cms_init'
            _init, _root, pkg_ver = find_installed_pkg(pkg, debug)
            if  pkg == 'external/py2-matplotlib':
                matplotlib_ver = pkg_ver
            if  not pkg_ver:
                continue
            if  pkg == 'lcg/root':
                func = 'link_root'
            msg += '%s "%s" "%s"\n' % (func, pkg, pkg_ver)
        msg += '# Recreate release dirs in order to add them to PATHs\n'
        msg += 'mkdir -p $CMSSH_ROOT/install/lib/release_lib\n'
        msg += 'mkdir -p $CMSSH_ROOT/install/lib/release_root/lib\n'
        msg += 'mkdir -p $CMSSH_ROOT/install/lib/release_external/lib\n'
        msg += 'export DYLD_LIBRARY_PATH=$CMSSH_ROOT/globus/lib:$CMSSH_ROOT/glite/lib:$CMSSH_ROOT/install/lib\n'
        msg += 'export DYLD_LIBRARY_PATH=$CMSSH_ROOT/install/lib/release_lib:$CMSSH_ROOT/install/lib/release_external:$CMSSH_ROOT/install/lib/release_external/lib:$DYLD_LIBRARY_PATH\n'
        msg += 'export LD_LIBRARY_PATH=$CMSSH_ROOT/globus/lib:$CMSSH_ROOT/glite/lib:$CMSSH_ROOT/install/lib:$LD_LIBRARY_PATH\n'
        msg += 'export LD_LIBRARY_PATH=$CMSSH_ROOT/install/lib/release_lib:$CMSSH_ROOT/install/lib/release_external/lib:$CMSSH_ROOT/install/lib/release_root/lib:$LD_LIBRARY_PATH\n'
        if  parch == 'x86_64':
            msg += 'export LD_LIBRARY_PATH=$CMSSH_ROOT/globus/lib64:$CMSSH_ROOT/glite/lib64:$CMSSH_ROOT/install/lib64:$LD_LIBRARY_PATH\n'
        msg += 'export PATH=$VO_CMS_SW_DIR/bin:$CMSSH_ROOT/install/bin:$PATH\n'
        msg += 'export PATH=$PATH:$CMSSH_ROOT/globus/bin\n'
        msg += 'export PATH=$PATH:$CMSSH_ROOT/glite/bin\n'
        msg += 'export PATH=$PATH:$CMSSH_ROOT/srmclient2/bin\n'
        msg += 'export PATH=$PATH:$CMSSH_ROOT/bin\n'
        msg += 'export PATH=$PATH:$CMSSH_ROOT/lcg/bin\n'
        msg += 'export PATH=$PATH:$CMSSH_ROOT/CRABClient/bin\n'
        msg += 'set_cmd "lcg"\n'
        msg += 'set_cmd "srm"\n'
        msg += 'export PYTHONPATH=$CMSSH_ROOT/cmssh/src:$PYTHONPATH\n'
        msg += 'export PYTHONPATH=$CMSSH_ROOT/install/lib/release_root/lib:$PYTHONPATH\n'
        msg += 'export PYTHONPATH=$PYTHONPATH:$CMSSH_ROOT\n'
        msg += 'export PYTHONPATH=$PYTHONPATH:$CMSSH_ROOT/CRABClient/src/python\n'
        msg += 'export PYTHONPATH=$PYTHONPATH:$CMSSH_ROOT/WMCore/src/python\n'
        msg += 'export PYTHONPATH=$CMSSH_ROOT/install/lib/python%s/site-packages:$PYTHONPATH\n' % pver
        msg += 'export DBS_INSTANCE=cms_dbs_prod_global\n'
        msg += 'export LCG_GFAL_INFOSYS=lcg-bdii.cern.ch:2170\n'
        msg += 'export VOMS_USERCONF=$CMSSH_ROOT/glite/etc/vomses\n'
        msg += 'export VOMS_LOCATION=$CMSSH_ROOT/glite\n'
        msg += 'export MYPROXY_SERVER=myproxy.cern.ch\n'
        msg += 'export X509_CERT_DIR=$CMSSH_ROOT/certificates\n'
        msg += 'export GLOBUS_ERROR_VERBOSE=true\n'
        msg += 'export GLOBUS_OPTIONS=-Xmx512M\n'
        msg += 'export GLOBUS_TCP_PORT_RANGE=34000,35000\n'
        msg += 'export GLOBUS_PATH=$CMSSH_ROOT/globus\n'
        msg += 'export GLOBUS_LOCATION=$CMSSH_ROOT/globus\n'
        msg += 'export VOMS_PROXY_INFO_DONT_VERIFY_AC=anything_you_want\n'
        msg += 'export MATPLOTLIBRC=$CMSSH_ROOT/cmssh/src/config\n'
        msg += 'export CORAL_DIR=$CMSSH_ROOT/CMSSW/$SCRAM_ARCH/cms/coral/%s\n' % coral_ver
        fname = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt'
        if  os.path.isfile(fname):
            msg += 'export CMS_JSON=%s\n' % fname
        else:
            msg += 'export CMS_JSON=None\n'
        msg += 'coral_init\n'
        cmssw = opts.cmssw if opts.cmssw else ''
        msg += 'export CMSSH_CMSSW=%s\n' % cmssw
        msg += 'echo " DONE"'
        if  debug:
            print "+++ write setup.sh"
        setup.write(msg)
    os.chmod('setup.sh', 0755)
    if  not os.path.islink(os.path.join(path, 'CMSSW')):
        sconfig_dir = os.path.join(path, 'CMSSW/SITECONF/local/JobConfig')
        try:
            os.makedirs(sconfig_dir)
        except:
            pass
        with open(os.path.join(sconfig_dir, 'site-local-config.xml'), 'w') as sconfig:
            sconfig.write(siteconfig() + '\n')

    vomses = os.path.join(path, 'glite')
    print "Create vomses area"
    vdir = os.path.join(vomses, 'etc')
    try:
        os.makedirs(vdir)
    except:
        pass
    fname = os.path.join(vdir, 'vomses')
    if  not os.path.isfile(fname):
        with open(fname, 'w') as fds:
            msg = '"cms" "voms.fnal.gov" "15015" "/DC=org/DC=doegrids/OU=Services/CN=http/voms.fnal.gov" "cms"'
            fds.write(msg + '\n')
            msg = '"cms" "voms.cern.ch" "15002" "/DC=ch/DC=cern/OU=computers/CN=voms.cern.ch" "cms"'
            fds.write(msg + '\n')
            msg = '"cms" "lcg-voms.cern.ch" "15002" "/DC=ch/DC=cern/OU=computers/CN=lcg-voms.cern.ch" "cms"'
            fds.write(msg + '\n')
        os.chmod(fname, stat.S_IRUSR | stat.S_IROTH | stat.S_IRGRP)

    print "Create cmssh"
    try:
        os.makedirs(os.path.join(path, 'bin'))
    except:
        pass
    with open(os.path.join(path, 'bin/cmssh'), 'w') as cmssh:
        msg  = '#!/bin/bash\n'
        msg += """
if  [ $# == 1 ]; then
    list="help -help --help -h"
    if [[ $list =~ $1 ]]; then
        echo "Usage: $0 <notebook>"
        echo "      notebook - start cmssh in notebook mode"
        echo "                 (will start cmssh session in a browser)"
        exit;
    fi
fi\n"""
        msg += 'echo "Welcome to cmssh, %s@%s"\n' % (cmssh_ver, cmssh_ts)
        msg += 'source %s/setup.sh\n' % path
        if  opts.multi_user:
            msg += 'ipdir="/tmp/$USER/.ipython"\nmkdir -p $ipdir\n'
        else:
            msg += 'ipdir="%s/.ipython"\nmkdir -p $ipdir\n' % path
        msg += """
soft_dir=%(path)s
if [ ! -d $ipdir/extensions ]; then
    mkdir -p $ipdir/extensions
fi
if [ ! -d $ipdir/profile_cmssh ]; then
    mkdir -p $ipdir/profile_cmssh
fi
if [ -f $ipdir/extensions/cmssh_extension.py ]; then
    /bin/rm -f /tmp/$USER/.ipython/extensions/cmssh_extension.py
fi
if [ ! -f $ipdir/extensions/cmssh_extension.py ]; then
    cp $soft_dir/cmssh/src/config/cmssh_extension.py $ipdir/extensions/
fi
if [ ! -f $ipdir/profile_cmssh/ipython_config.py ]; then
    cp $soft_dir/cmssh/src/config/ipython_config.py $ipdir/profile_cmssh/
fi
if [ ! -f $HOME/.globus/userkey.pem ]; then
    echo "You don't have $HOME/.globus/userkey.pem on this system"
    echo "Please install it to proceed"
    exit -1
fi
if [ ! -f $HOME/.globus/usercert.pem ]; then
    echo "You don't have $HOME/.globus/usercert.pem on this system"
    echo "Please install it to proceed"
    exit -1
fi
export IPYTHONDIR=$ipdir
pylab=" --pylab=auto"
osname=`uname -s`
if  [ "$osname" == "Darwin" ]; then
    cms_osx_driver=$CMSSH_ROOT/CMSSW/$SCRAM_ARCH/external/py2-matplotlib/%(mver)s/lib/python2.7/site-packages/matplotlib/backends_macosx.so
    osx_driver=$CMSSH_ROOT/install/lib/python%(pver)s/site-packages/matplotlib/backends/_macosx.so
    if  [ -f "$osx_driver" ] || [ -f "$cms_osx_driver" ]; then
        pylab=" --pylab=osx"
    fi
fi
notebook="--no-banner"
if  [ $# == 1 ]; then
    if [ $1 == "notebook" ]; then
        notebook="notebook"
        pylab="--pylab=inline"
        export CMSSH_NOTEBOOK=1
    fi
fi
opts="$notebook $pylab"
ipython $opts --ipython-dir=$ipdir --profile=cmssh
""" % {'path':path, 'pver':pver, 'mver': matplotlib_ver}
        cmssh.write(msg)
    os.chmod('bin/cmssh', 0755)

    # remove soflinks in coral, since matplotlib fails to load
    if  platform == 'Darwin' and osx_ver() == '10.6':
        cdir = coral_root + '/external/%s/lib' % get_scram_arch()
        for item in os.listdir(cdir):
            fname = os.path.join(cdir, item)
            if  os.path.islink(fname):
                os.remove(fname)

    print "Make links"
    xrdcp = os.path.join(path, 'install/bin/xrdcp')
    if  not os.path.islink(xrdcp):
        _xrootd_init, xrootd_root, _xrootd_ver = find_installed_pkg('external/xrootd', debug)
        os.symlink(os.path.join(xrootd_root, 'bin/xrdcp'), xrdcp)

    print "Clean-up soft area"
    os.chdir(path)
    subprocess.call("rm *.tar.gz", shell=True)
    if  not os.path.isdir('logs'):
        os.makedirs('logs')
    subprocess.call("mv *.log logs", shell=True)

    print "Congratulations, cmssh is available at %s/bin/cmssh" % path

if __name__ == '__main__':
    main()
