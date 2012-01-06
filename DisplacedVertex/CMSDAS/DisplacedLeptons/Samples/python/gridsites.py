# list of sites obtained with
# dbs search --query="find site"

gridsites_good=["heplnx204.pp.rl.ac.uk","gfe02.hep.ph.ic.ac.uk","dcache-se-cms.desy.de",
                "dc2-grid-64.brunel.ac.uk","gfe02.hep.ph.ic.ac.uk","grid-ce.physik.rwth-aachen.de"]
gridsites_bad=["ucsd.edu", # low throughput
               "ufl.edu",  # low throughput
               "indiacms", # all jobs aborted
               "metu.edu.tr", # all jobs aborted
               "srm-cms.gridpp.rl.ac.uk", # not allowed to use (Tier1)
               "caf.cern.ch", # not allowed
               "cmssrm-fzk.gridka.de", # not allowed: T1_DE_KIT
               "cmssrm.fnal.gov", # not allowed
               "lcg38.sinp.msu.ru", # exit code 8001
               "grid.sinica.edu.tw", # jobs stay scheduled forever
               "srm.unl.edu",
               "cmsdca2.fnal.gov" # not allowed to use (FNAL tier3)
               ]

# full list of sites below. most of which I have no opinion about (yet)
# =====================================================================
#
# uscms1-se.fltech-grid3.fit.edu
# tigressdata.princeton.edu
# t3se01.psi.ch
# t2se01.physics.ox.ac.uk
# t2-srm-02.lnl.infn.it
# svr018.gla.scotgrid.ac.uk
# storm.ifca.es
# storm-se-01.ba.infn.it
# storm-fe-cms.cr.cnaf.infn.it
# storage01.lcg.cscs.ch
# srmcms.pic.es
# srm2.grid.sinica.edu.tw
# srm01.ncg.ingrid.pt
# srm01.lip.pt
# srm.unl.edu
# srm.ucr.edu
# srm.test3.edu
# srm.princeton.edu
# srm.minnesota.edu
# srm.ihepa.ufl.edu
# srm.ihep.ac.cn
# srm.glite.ecdf.ed.ac.uk
# srm.ciemat.es
# srm-dcache.rcac.purdue.edu
# srm-cms.gridpp.rl.ac.uk
# srm-cms.cern.ch
# sigmorgh.hpcc.ttu.edu
# se3.itep.ru
# se2.ppgrid1.rhul.ac.uk
# se1.accre.vanderbilt.edu
# se03.esc.qmul.ac.uk
# se01.indiacms.res.in
# se01.cmsaf.mit.edu
# se0003.m45.ihep.su
# se.tier3.ucdavis.edu
# se.polgrid.pl
# se.grid.kiae.ru
# se-dcache.hepgrid.uerj.br
# sbgse1.in2p3.fr
# ruhex-osgce.rutgers.edu
# polgrid4.in2p3.fr
# phedex.geol.uniovi.es
# pcncp22.ncp.edu.pk
# pc138.physics.uoi.gr
# osg-se.sprace.org.br
# osg-se.cac.cornell.edu
# osg-hep.phys.virginia.edu
# node12.datagrid.cea.fr
# meson.fis.cinvestav.mx
# maite.iihe.ac.be
# madhatter.csc.fi
# lyogrid06.in2p3.fr
# lcgsedc01.jinr.ru
# lcgse02.phy.bris.ac.uk
# io.hep.kbfi.ee
# ingrid-se02.cism.ucl.ac.be
# hepse01.colorado.edu
# heplnx204.pp.rl.ac.uk
# hephyse.oeaw.ac.at
# hepcms-0.umd.edu
# hep.pha.jhu.edu
# gw-3.ccc.ucl.ac.uk
# grse001.inr.troitsk.ru
# gridse3.pg.infn.it
# grid143.kfki.hu
# grid128.sinp.msu.ru
# grid02.phy.ncu.edu.tw
# grid-srm.physik.rwth-aachen.de
# gfe02.hep.ph.ic.ac.uk
# ff-se.unl.edu
# f-dpm001.grid.sinica.edu.tw
# eymir.grid.metu.edu.tr
# dcache-se-cms.desy.de
# dc2-grid-64.brunel.ac.uk
# cmssrm.hep.wisc.edu
# cmssrm.fnal.gov
# cmssrm-fzk.gridka.de
# cmsrm-se01.roma1.infn.it
# cmsdcache.pi.infn.it
# cmsdca2.fnal.gov
# cms-xen19.fnal.gov
# cms-se0.kipt.kharkov.ua
# cms-0.mps.ohio-state.edu
# cluster142.knu.ac.kr
# cluster.pnpi.nw.ru
# cit-se.ultralight.org
# ccsrmt2.in2p3.fr
# ccsrm.in2p3.fr
# caf.cern.ch
# bsrm-1.t2.ucsd.edu
# bonner-cms3.rice.edu
