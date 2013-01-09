class check_HN_name:
    def init(self): 
        pass
    
    def stdaloneCheck(self):  
 
        import urllib
        from commands import getstatusoutput
        print 'start standalone check ...\n'
        status, dn = getstatusoutput('voms-proxy-info -identity')
        if status == 0:
           print "my DN is: %s \n"%dn
        dn = dn.split('\n')[-1]
        dn  = urllib.urlencode({'dn':dn})
        print 'Using urlencoded DN: \n\t %s '%dn
        f = urllib.urlopen("https://cmsweb.cern.ch/sitedb/json/index/dnUserName?%s" % dn)
        print 'my HN user name is: %s \n'%str(f.read())
        f.close()
        print '\nend check.....................'

    def crabCheck(self):
        from CrabLogger import CrabLogger
        from WorkSpace import *
        import tempfile, urllib, os 
      
        dname = tempfile.mkdtemp( "", "crab_", '/tmp' )
        os.system("mkdir %s/log"%dname )
        os.system("touch %s/crab.log"%dname )
        
        cfg_params={'USER.logdir' : dname }
        common.work_space = WorkSpace(dname, cfg_params)
        args = string.join(sys.argv,' ')
        common.debugLevel = 0
        common.logger = CrabLogger(args)
        
        from crab_util import getDN,gethnUserNameFromSiteDB
        print 'start using CRAB utils ...\n'
        print "my DN is: %s \n"%getDN()
        try:
            print 'my HN user name is: %s \n'%gethnUserNameFromSiteDB()
        except:
            print '\nWARNING native crab_utils failed! ' 
            dn=urllib.urlencode({'dn':getDN()})
            print 'trying now using urlencoded DN: \n\t %s '%dn
            status,hnName = self.gethnName_urlenc(dn)
            if status == 1: 
                print '\nWARNING: failed also using urlencoded DN '
            else: 
                print 'my HN user name is: %s \n'%name
                print 'problems with crab_utils'   
        print '\nend check.....................'
        
        os.system("rm -rf %s"%dname )
         
    def gethnName_urlenc(self,dn):
        from WMCore.Services.SiteDB.SiteDB import SiteDBJSON
        hnUserName = None
        userdn = dn
        mySiteDB = SiteDBJSON()
        status = 0 
        try:
            hnUserName = mySiteDB.dnUserName(dn=userdn)
        except:
            status = 1 
        return status,hnUserName


if __name__ == '__main__' :
    import sys
    args = sys.argv[1:]
    check = check_HN_name() 
    if 'crab' in args:
        check.crabCheck()  
    else:
        check.stdaloneCheck()  

