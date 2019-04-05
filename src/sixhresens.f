      program sixhresens

c$$$ to compile, use Makefile (need netcdf.mod, wrf_tools.mod in dir, compile
c$$$ module_netcdf.f and module_wrf_tools.f if the .mod files don't exist)
c$$$ Also, remove enscalc before compiling with Makefile.  Lastly, make sure to
c$$$ copy a single ens member over to all of the outfiles (see below)

c$$$ This program reads data from fens files and anal files and
c$$$ calculates ensemble senstivity values for six hour response functions

c$$$ Originally by Brian Ancell; adapted for six-hour response functions
c$$  by Austin Coleman.
c$$  11/20/2018

      USE NETCDF
      USE WRF_TOOLS

      real rp,rpmean,latapp,lonapp,igg,jgg
      type(proj_info) :: proj_in
      character*67 filectw,filec,modfile
      character*67 outfile1
      integer fidc,idvar,cntall,ii
      real ir,jr,missing2
      integer idex,jdex

      real, pointer::dat(:,:)
      real, pointer::dat2(:,:)
      real, pointer::ucross(:,:,:)
      real, pointer::vcross(:,:,:)
      real, pointer::ucrossmean(:,:,:)
      real, pointer::vcrossmean(:,:,:)
      integer mix,mjx,mkx,permissr,permissw,rnum,iunitobs
      integer iunitwrite1,iunitread,timesens,timeresponse
      integer i,j,k,f
      real q2max, q2maxii, q2maxjj, u10maxii, u10maxjj, u10max
      real, pointer::gph300mean(:,:)
      real, pointer::gph500mean(:,:)
      real, pointer::gph700mean(:,:)
      real, pointer::gph850mean(:,:)
      real, pointer::gph925mean(:,:)
      real, pointer::u300mean(:,:)
      real, pointer::u500mean(:,:)
      real, pointer::u700mean(:,:)
      real, pointer::u850mean(:,:)
      real, pointer::u925mean(:,:)
      real, pointer::t300mean(:,:)
      real, pointer::t500mean(:,:)
      real, pointer::t700mean(:,:)
      real, pointer::t850mean(:,:)
      real, pointer::t925mean(:,:)
      real, pointer::v300mean(:,:)
      real, pointer::v500mean(:,:)
      real, pointer::v700mean(:,:)
      real, pointer::v850mean(:,:)
      real, pointer::v925mean(:,:)
      real, pointer::gph300(:,:)
      real, pointer::gph500(:,:)
      real, pointer::gph700(:,:)
      real, pointer::gph850(:,:)
      real, pointer::gph925(:,:)
      real, pointer::u300(:,:)
      real, pointer::u500(:,:)
      real, pointer::u700(:,:)
      real, pointer::u850(:,:)
      real, pointer::u925(:,:)
      real, pointer::t300(:,:)
      real, pointer::t500(:,:)
      real, pointer::t700(:,:)
      real, pointer::t850(:,:)
      real, pointer::t925(:,:)
      real, pointer::v300(:,:)
      real, pointer::v500(:,:)
      real, pointer::v700(:,:)
      real, pointer::v850(:,:)
      real, pointer::v925(:,:)
      real, pointer::u10(:,:)
      real, pointer::v10(:,:)
      real, pointer::dbzmean(:,:,:)
      real, pointer::wspmean(:,:)
      real, pointer::uphelmean(:,:)
      real, pointer::dbzresponse(:,:,:)
      real, pointer::wspresponse(:,:)
      real, pointer::uphelresponse(:,:)
      real, pointer::dbz(:,:,:)
      real, pointer::wsp(:,:)
      real, pointer::uphel(:,:)
      real, pointer::qvapormean(:,:,:)
      real, pointer::totpresmean(:,:,:)
      real, pointer::tmean(:,:,:)
      real, pointer::phmean(:,:,:)
      real, pointer::phbmean(:,:,:)
      real, pointer::thetamean(:,:,:)
      real, pointer::tempmean(:,:,:)
      real, pointer::rhodrymean(:,:,:)
      real, pointer::uumean(:,:,:)
      real, pointer::vvmean(:,:,:)
      real, pointer::gpttotmean(:,:,:)
      real, pointer::gpheightmean(:,:,:)
      real, pointer::gpheightmeanhalf(:,:,:)
      real, pointer::pmean(:,:,:)
      real, pointer::pbmean(:,:,:)
      real, pointer::surfptotalmean(:,:)
      real, pointer::mumean(:,:)
      real, pointer::mubmean(:,:)
      real, pointer::psfcmean(:,:)
      real, pointer::t2mean(:,:)
      real, pointer::q2mean(:,:)
      real, pointer::slpmean(:,:)
      real, pointer::u10mean(:,:)
      real, pointer::v10mean(:,:)

      real, pointer::qvaporgfs(:,:,:)
      real, pointer::totpresgfs(:,:,:)
      real, pointer::tgfs(:,:,:)
      real, pointer::phgfs(:,:,:)
      real, pointer::phbgfs(:,:,:)
      real, pointer::thetagfs(:,:,:)
      real, pointer::tempgfs(:,:,:)
      real, pointer::rhodrygfs(:,:,:)
      real, pointer::uugfs(:,:,:)
      real, pointer::vvgfs(:,:,:)
      real, pointer::uucrossgfs(:,:,:)
      real, pointer::vvcrossgfs(:,:,:)
      real, pointer::gpttotgfs(:,:,:)
      real, pointer::gpheightgfs(:,:,:)
      real, pointer::gpheighthalfgfs(:,:,:)
      real, pointer::surfptotalgfs(:,:)
      real, pointer::mugfs(:,:)
      real, pointer::mubgfs(:,:)
      real, pointer::psfcgfs(:,:)
      real, pointer::t2gfs(:,:)
      real, pointer::q2gfs(:,:)
      real, pointer::slpgfs(:,:)
      real, pointer::gph300gfs(:,:)
      real, pointer::gph500gfs(:,:)
      real, pointer::gph700gfs(:,:)
      real, pointer::gph850gfs(:,:)
      real, pointer::t850gfs(:,:)
      real, pointer::evec1(:)

      real, pointer::qvaporresponse(:,:,:)
      real, pointer::totpresresponse(:,:,:)
      real, pointer::tresponse(:,:,:)
      real, pointer::phresponse(:,:,:)
      real, pointer::phbresponse(:,:,:)
      real, pointer::thetaresponse(:,:,:)
      real, pointer::tempresponse(:,:,:)
      real, pointer::rhodryresponse(:,:,:)
      real, pointer::gpttotresponse(:,:,:)
      real, pointer::gpheightresponse(:,:,:)
      real, pointer::gpheighthalfresponse(:,:,:)
      real, pointer::surfptotalresponse(:,:)
      real, pointer::muresponse(:,:)
      real, pointer::mubresponse(:,:)
      real, pointer::psfcresponse(:,:)
      real, pointer::t2response(:,:)
      real, pointer::q2response(:,:)
      real, pointer::slpresponse(:,:)
      real, pointer::u10response(:,:)
      real, pointer::v10response(:,:)
      real, pointer::t850response(:,:)
      real, pointer::searainresponsef(:,:)
      real, pointer::searainresponsei(:,:)
      real, pointer::uuvar(:,:,:)
      real, pointer::vvvar(:,:,:)
      real, pointer::uucrossvar(:,:,:)
      real, pointer::vvcrossvar(:,:,:)
      real, pointer::uucross(:,:,:)
      real, pointer::vvcross(:,:,:)
      real, pointer::uucrossmean(:,:,:)
      real, pointer::vvcrossmean(:,:,:)
      real, pointer::uucrosscovar(:,:,:)
      real, pointer::vvcrosscovar(:,:,:)
      real, pointer::tempvar(:,:,:)
      real, pointer::totpresvar(:,:,:)
      real, pointer::slpvar(:,:)
      real, pointer::psfcvar(:,:)

      real, pointer::u300var(:,:)
      real, pointer::u500var(:,:)
      real, pointer::u700var(:,:)
      real, pointer::u850var(:,:)
      real, pointer::u925var(:,:)
      real, pointer::t300var(:,:)
      real, pointer::t500var(:,:)
      real, pointer::t700var(:,:)
      real, pointer::t850var(:,:)
      real, pointer::t925var(:,:)
      real, pointer::v300var(:,:)
      real, pointer::v500var(:,:)
      real, pointer::v700var(:,:)
      real, pointer::v850var(:,:)
      real, pointer::v925var(:,:)
      real, pointer::gph300var(:,:)
      real, pointer::gph500var(:,:)
      real, pointer::gph700var(:,:)
      real, pointer::gph850var(:,:)
      real, pointer::gph925var(:,:)
      real, pointer::q850var(:,:)
      real, pointer::q850mean(:,:)
      real, pointer::q850(:,:)
      real, pointer::t2var(:,:)
      real, pointer::td2var(:,:)
      real, pointer::td2mean(:,:)
      real, pointer::td2(:,:)
      real, pointer::q2var(:,:)
      real, pointer::u10var(:,:)
      real, pointer::v10var(:,:)

      real, pointer::gph300r1covar(:,:,:)
      real, pointer::gph500r1covar(:,:,:)
      real, pointer::gph700r1covar(:,:,:)
      real, pointer::gph850r1covar(:,:,:)
      real, pointer::gph925r1covar(:,:,:)
      real, pointer::gph300targ(:,:,:)
      real, pointer::gph500targ(:,:,:)
      real, pointer::gph700targ(:,:,:)
      real, pointer::gph850targ(:,:,:)
      real, pointer::gph925targ(:,:,:)

      real, pointer::u300r1covar(:,:,:)
      real, pointer::u500r1covar(:,:,:)
      real, pointer::u700r1covar(:,:,:)
      real, pointer::u850r1covar(:,:,:)
      real, pointer::u925r1covar(:,:,:)
      real, pointer::v300r1covar(:,:,:)
      real, pointer::v500r1covar(:,:,:)
      real, pointer::v700r1covar(:,:,:)
      real, pointer::v850r1covar(:,:,:)
      real, pointer::v925r1covar(:,:,:)
      real, pointer::t300r1covar(:,:,:)
      real, pointer::t500r1covar(:,:,:)
      real, pointer::t700r1covar(:,:,:)
      real, pointer::t850r1covar(:,:,:)
      real, pointer::t925r1covar(:,:,:)

      real, pointer::u300targ(:,:,:)
      real, pointer::u500targ(:,:,:)
      real, pointer::u700targ(:,:,:)
      real, pointer::u850targ(:,:,:)
      real, pointer::u925targ(:,:,:)
      real, pointer::v300targ(:,:,:)
      real, pointer::v500targ(:,:,:)
      real, pointer::v700targ(:,:,:)
      real, pointer::v850targ(:,:,:)
      real, pointer::v925targ(:,:,:)
      real, pointer::t300targ(:,:,:)
      real, pointer::t500targ(:,:,:)
      real, pointer::t700targ(:,:,:)
      real, pointer::t850targ(:,:,:)
      real, pointer::t925targ(:,:,:)

      real, pointer::q2r1covar(:,:,:)
      real, pointer::t2r1covar(:,:,:)
      real, pointer::td2r1covar(:,:,:)
      real, pointer::q850r1covar(:,:,:)
      real, pointer::u10r1covar(:,:,:)
      real, pointer::v10r1covar(:,:,:)

      real, pointer::slpr1covar(:,:,:)

      real, pointer::q2targ(:,:,:)
      real, pointer::t2targ(:,:,:)
      real, pointer::td2targ(:,:,:)
      real, pointer::q850targ(:,:,:)
      real, pointer::u10targ(:,:,:)
      real, pointer::v10targ(:,:,:)

      real, pointer::slptarg(:,:,:)

      real, pointer::qvapor(:,:,:)
      real, pointer::totpres(:,:,:)
      real, pointer::t(:,:,:)
      real, pointer::ph(:,:,:)
      real, pointer::phb(:,:,:)
      real, pointer::theta(:,:,:)
      real, pointer::temp(:,:,:)
      real, pointer::rhodry(:,:,:)
      real, pointer::uu(:,:,:)
      real, pointer::vv(:,:,:)
      real, pointer::gpttot(:,:,:)
      real, pointer::gpheight(:,:,:)
      real, pointer::gpheighthalf(:,:,:)
      real, pointer::p(:,:,:)
      real, pointer::pb(:,:,:)
      real, pointer::surfptotal(:,:)
      real, pointer::mu(:,:)
      real, pointer::mub(:,:)
      real, pointer::psfc(:,:)
      real, pointer::t2(:,:)
      real, pointer::q2(:,:)
      real, pointer::hgt(:,:)
      real, pointer::slp(:,:)
      real, pointer::pcptotmean(:,:)
      real, pointer::pcptotresponse(:,:)
      real, pointer::pcptot(:,:)
      real, pointer::smat(:,:,:)
      real, pointer::tmat(:,:,:)
      real, pointer::mmat(:,:,:)

      integer, pointer::uggph850(:,:)
      integer, pointer::uggph925(:,:)
      integer, pointer::uggph700(:,:)
      integer iunit1
      real, pointer::znu(:)
      real, pointer::znw(:)

      character*67, pointer::infileSENS(:)
      character*67, pointer::infileR(:,:)
      character*67, pointer::infile(:)
      real basethet,levee,rdgas,rvgas,preffer
      real kap,prefcont,kapdiv
      integer m,ensnum,ensnumax
      integer numpoints
      character*20 infilemeanSENS,infilemeanR, infilerefvar
      character*20 infilerefvard2,outfile2
      integer iunitmean, iunitref, rcode,numtimes,timer, iunit2
      real grav,obdiff,obspd,pmin
      integer cyccenti,cyccentj,numcyc,numgrad
      real grad1mean,grad2mean,grad3mean,grad4mean
      real grad1,grad2,grad3,grad4,rumean,rvmean,ruu,rvv
      integer is,ie,js,je
      integer RR,numresp,uhthresh,dbzthresh
      character*9 filer
      integer iddbzavg,iddbzmax,iduhavg,iduhmax,idpcp,idwspd,iduhcov,mem
      integer iddbzcov

      real dbzavg
      real windavg
      real pcp
      real dbzmax
      real uhavg
      real uhmax
      real uhcov
      real dbzcov
      real LONis,LONie,LATjs,LATje

      real dbzavgMEM
      real windavgMEM
      real pcpMEM
      real dbzmaxMEM
      real uhavgMEM
      real uhmaxMEM
      real uhcovMEM
      real dbzcovMEM
      real gobvar300,gobvar500,gobvar700,gobvar850
      real gobvar925
      real tobvar300,tobvar500,tobvar700,tobvar850
      real tobvar925,qobvar850,w10obvar
      real wobvar300,wobvar500,wobvar700,wobvar850
      real wobvar925,t2obvar,q2obvar,slpobvar,td2obvar
      integer rrr
      character*2 senstime
      character*2 rtime
      character*2 time
      logical scatter, exists
      integer sensvar
      integer rix,rjx,rkx
      integer fidr, npts
      real fmiss

      real, pointer :: dbzmaxvec(:)
      real, pointer :: dbzavgvec(:)
      real, pointer :: uhavgvec(:)
      real, pointer :: uhmaxvec(:)
      real, pointer :: pcpvec(:)
      real, pointer :: windavgvec(:)
      real, pointer :: uhcovvec(:)
      real, pointer :: dbzcovvec(:)

      real, pointer::qvapornext(:,:,:)
      real, pointer::dbznext(:,:,:)
      real, pointer::totpresnext(:,:,:)
      real, pointer::uhnext(:,:)
      real, pointer::tnext(:,:,:)
      real, pointer::phnext(:,:,:)
      real, pointer::phbnext(:,:,:)
      real, pointer::thetanext(:,:,:)
      real, pointer::tempnext(:,:,:)
      real, pointer::rhodrynext(:,:,:)
      real, pointer::uunext(:,:,:)
      real, pointer::vvnext(:,:,:)
      real, pointer::wwnext(:,:,:)
      real, pointer::gpttotnext(:,:,:)
      real, pointer::pnext(:,:,:)
      real, pointer::phydnext(:,:,:)
      real, pointer::pbnext(:,:,:)
      real, pointer::qrainnext(:,:,:)
      real, pointer::qcloudnext(:,:,:)
      real, pointer::qgraupnext(:,:,:)
      real, pointer::qicenext(:,:,:)
      real, pointer::qsnownext(:,:,:)
      real, pointer::surfptotalnext(:,:)
      real, pointer::munext(:,:)
      real, pointer::mubnext(:,:)
      real, pointer::sstnext(:,:)
      real, pointer::rainncnext(:,:)
      real, pointer::snowncnext(:,:)
      real, pointer::hailncnext(:,:)
      real, pointer::graupelncnext(:,:)
      real, pointer::rainnumnext(:,:,:)
      real, pointer::icenumnext(:,:,:)
      real, pointer::q2next(:,:)
      real, pointer::t2next(:,:)
      real, pointer::th2next(:,:)
      real, pointer::psfcnext(:,:)
      real, pointer::u10next(:,:)
      real, pointer::v10next(:,:)
      real, pointer::rdnwnext(:)
      real, pointer::rdnnext(:)
      real, pointer::dnwnext(:)
      real, pointer::dnnext(:)
      real, pointer::znunext(:)
      real, pointer::znwnext(:)
      real, pointer::wspdnext(:,:)
      real uhmaxmean
      real dbzmaxmean
      real uhcovmean
      real dbzcovmean

      integer iunitwritemean,cc
      integer, pointer::levels(:)
      character*20 outfilemean
      character*2 fhr, ensmem

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C$$$ End of Declarations $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

C$$$$ read in namelist values for R calculations

      read*,ensnum
      read*,timesens
      read*,timer
      read*,rrr
      read*,LONis
      read*,LONie
      read*,LATjs
      read*,LATje

c$$ Optional add on - save member response vectors to netCDF
c$$  by setting scatter to true

      read*,scatter

c$$ UH Threshold (only for UH Coverage)

      read*,uhthresh

c$$ DBZ Threshold (only for DBZ Coverage)

      read*,dbzthresh

c$$$ Initialize constants

      infilemeanSENS='SENSmean.out'
      infilemeanR='Rmean.out'
      infilerefvar='wrfoutREF'
      infilerefvard2='wrfoutREFd2'
      outfile1='wrfout.sens'

      ensnumax=42
      iunitread=10
      iunitmean=44
      iunitref = 32
      iunit1=23
      iunit2 = 24
      basethet=300.
      rdgas=287.
      rvgas=461.6
      preffer=100000.
      kap=2./7.
      grav=9.81
      prefcont=preffer**(-(kap))
      kapdiv=1./(1.-kap)
      permissw=3
      permissr=2
      missing2=-9999.0
      fmiss=missing2
      fidc=93
      timeresponse=timer
      numtimes=48

C$$$ Convert lat/lon to i,j

      call set_domain_proj(infilerefvard2, proj_in)

      call latlon_to_ij(proj_in, LATjs,
     &       LONis , ir, jr)
      is=floor(ir)
      js=floor(jr)
      call latlon_to_ij(proj_in, LATje,
     &       LONie , ir, jr)
      ie=floor(ir)
      je=floor(jr)

      print*, is,ie,js,je

c$$$  variable arrays

      allocate(infileSENS(ensnumax))
      allocate(infileR(6,ensnumax))

      write(senstime, '(i2)') timesens
      write(rtime, '(i2)') timer

      infileSENS = (/
     &     'mem1/SENS1_' // trim(adjustl(senstime)) // '.out',
     &     'mem2/SENS2_' // trim(adjustl(senstime)) // '.out',
     &     'mem3/SENS3_' // trim(adjustl(senstime)) // '.out',
     &     'mem4/SENS4_' // trim(adjustl(senstime)) // '.out',
     &     'mem5/SENS5_' // trim(adjustl(senstime)) // '.out',
     &     'mem6/SENS6_' // trim(adjustl(senstime)) // '.out',
     &     'mem7/SENS7_' // trim(adjustl(senstime)) // '.out',
     &     'mem8/SENS8_' // trim(adjustl(senstime)) // '.out',
     &     'mem9/SENS9_' // trim(adjustl(senstime)) // '.out',
     &     'mem10/SENS10_' // trim(adjustl(senstime)) // '.out',
     &     'mem11/SENS11_' // trim(adjustl(senstime)) // '.out',
     &     'mem12/SENS12_' // trim(adjustl(senstime)) // '.out',
     &     'mem13/SENS13_' // trim(adjustl(senstime)) // '.out',
     &     'mem14/SENS14_' // trim(adjustl(senstime)) // '.out',
     &     'mem15/SENS15_' // trim(adjustl(senstime)) // '.out',
     &     'mem16/SENS16_' // trim(adjustl(senstime)) // '.out',
     &     'mem17/SENS17_' // trim(adjustl(senstime)) // '.out',
     &     'mem18/SENS18_' // trim(adjustl(senstime)) // '.out',
     &     'mem19/SENS19_' // trim(adjustl(senstime)) // '.out',
     &     'mem20/SENS20_' // trim(adjustl(senstime)) // '.out',
     &     'mem21/SENS21_' // trim(adjustl(senstime)) // '.out',
     &     'mem22/SENS22_' // trim(adjustl(senstime)) // '.out',
     &     'mem23/SENS23_' // trim(adjustl(senstime)) // '.out',
     &     'mem24/SENS24_' // trim(adjustl(senstime)) // '.out',
     &     'mem25/SENS25_' // trim(adjustl(senstime)) // '.out',
     &     'mem26/SENS26_' // trim(adjustl(senstime)) // '.out',
     &     'mem27/SENS27_' // trim(adjustl(senstime)) // '.out',
     &     'mem28/SENS28_' // trim(adjustl(senstime)) // '.out',
     &     'mem29/SENS29_' // trim(adjustl(senstime)) // '.out',
     &     'mem30/SENS30_' // trim(adjustl(senstime)) // '.out',
     &     'mem31/SENS31_' // trim(adjustl(senstime)) // '.out',
     &     'mem32/SENS32_' // trim(adjustl(senstime)) // '.out',
     &     'mem33/SENS33_' // trim(adjustl(senstime)) // '.out',
     &     'mem34/SENS34_' // trim(adjustl(senstime)) // '.out',
     &     'mem35/SENS35_' // trim(adjustl(senstime)) // '.out',
     &     'mem36/SENS36_' // trim(adjustl(senstime)) // '.out',
     &     'mem37/SENS37_' // trim(adjustl(senstime)) // '.out',
     &     'mem38/SENS38_' // trim(adjustl(senstime)) // '.out',
     &     'mem39/SENS39_' // trim(adjustl(senstime)) // '.out',
     &     'mem40/SENS40_' // trim(adjustl(senstime)) // '.out',
     &     'mem41/SENS41_' // trim(adjustl(senstime)) // '.out',
     &     'mem42/SENS42_' // trim(adjustl(senstime)) // '.out' /)

      infileR(:,:) = ''
      print*, rtime
      f = timer
      do cc=timer-5,timer
          print*, cc
          do m=1,ensnum
              write(time, '(i2)') cc
              write(ensmem, '(i2)') m
c$$   get t back into i=1,2,3,4,5,6 etc. form
c$$    by adding six and subtracting final hr
c$$    print*, t-f+6
             infileR(cc-f+6,m) = ('mem' //trim(adjustl(ensmem)) //
     &        '/R' // trim(adjustl(ensmem)) // '_' //
     &        trim(adjustl(time)) // '.out')
              print*, infileR(cc-f+6, m)
          enddo
      enddo

c$$$  Get sens grid dims
      call open_file(infilerefvar,permissr,iunitref)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'WEST-EAST_GRID_DIMENSION', mix)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'SOUTH-NORTH_GRID_DIMENSION', mjx)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'BOTTOM-TOP_GRID_DIMENSION', mkx)
      allocate(znw(mkx))
      allocate(znu(mkx-1))
      call get_variable1d(iunitref,'ZNW',mkx,
     &     1,znw)
      call get_variable1d(iunitref,'ZNU',mkx-1,
     &     1,znu)
      call close_file(iunitref)

c$$$  Get response grid dims
      call open_file(infilerefvard2,permissr,iunitref)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'WEST-EAST_GRID_DIMENSION', rix)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'SOUTH-NORTH_GRID_DIMENSION', rjx)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'BOTTOM-TOP_GRID_DIMENSION', rkx)
      call close_file(iunitref)

      rkx = mkx

      allocate(evec1(mkx-1))
      allocate(smat(mix-1,mjx-1,mkx-1))
      allocate(tmat(mix-1,mjx-1,mkx-1))
      allocate(mmat(mix-1,mjx-1,mkx-1))

      allocate(uumean(mix,mjx-1,mkx-1))
      allocate(vvmean(mix-1,mjx,mkx-1))
      allocate(uucrossmean(mix-1,mjx-1,mkx-1))
      allocate(vvcrossmean(mix-1,mjx-1,mkx-1))
      allocate(phmean(mix-1,mjx-1,mkx))
      allocate(phbmean(mix-1,mjx-1,mkx))
      allocate(tmean(mix-1,mjx-1,mkx-1))
      allocate(mumean(mix-1,mjx-1))
      allocate(mubmean(mix-1,mjx-1))
      allocate(qvapormean(mix-1,mjx-1,mkx-1))
      allocate(totpresmean(mix-1,mjx-1,mkx-1))
      allocate(thetamean(mix-1,mjx-1,mkx-1))
      allocate(dbzmean(rix-1,rjx-1,rkx-1))
      allocate(dbzresponse(rix-1,rjx-1,rkx-1))
      allocate(dbz(rix-1,rjx-1,rkx-1))
      allocate(tempmean(mix-1,mjx-1,mkx-1))
      allocate(gpttotmean(mix-1,mjx-1,mkx))
      allocate(gpheightmean(mix-1,mjx-1,mkx))
      allocate(gpheightmeanhalf(mix-1,mjx-1,mkx-1))
      allocate(rhodrymean(mix-1,mjx-1,mkx-1))
      allocate(surfptotalmean(mix-1,mjx-1))
      allocate(psfcmean(mix-1,mjx-1))
      allocate(t2mean(mix-1,mjx-1))
      allocate(uphelmean(mix-1,mjx-1))
      allocate(wspmean(mix-1,mjx-1))
      allocate(uphelresponse(rix-1,rjx-1))
      allocate(wspresponse(rix-1,rjx-1))
      allocate(wsp(rix-1,rjx-1))
      allocate(uphel(rix-1,rjx-1))
      allocate(q2mean(mix-1,mjx-1))
      allocate(u10mean(mix-1,mjx-1))
      allocate(v10mean(mix-1,mjx-1))
      allocate(u10(mix-1,mjx-1))
      allocate(v10(mix-1,mjx-1))
      allocate(slpmean(mix-1,mjx-1))
      allocate(hgt(mix-1,mjx-1))
      allocate(pcptotmean(rix-1,rjx-1))
      allocate(pcptotresponse(rix-1,rjx-1))
      allocate(pcptot(rix-1,rjx-1))

      allocate(uu(mix,mjx-1,mkx-1))
      allocate(vv(mix-1,mjx,mkx-1))
      allocate(ph(mix-1,mjx-1,mkx))
      allocate(phb(mix-1,mjx-1,mkx))
      allocate(t(mix-1,mjx-1,mkx-1))
      allocate(mu(mix-1,mjx-1))
      allocate(mub(mix-1,mjx-1))
      allocate(qvapor(mix-1,mjx-1,mkx-1))
      allocate(totpres(mix-1,mjx-1,mkx-1))
      allocate(theta(mix-1,mjx-1,mkx-1))
      allocate(temp(mix-1,mjx-1,mkx-1))
      allocate(gpttot(mix-1,mjx-1,mkx))
      allocate(gpheight(mix-1,mjx-1,mkx))
      allocate(gpheighthalf(mix-1,mjx-1,mkx-1))
      allocate(rhodry(mix-1,mjx-1,mkx-1))
      allocate(surfptotal(mix-1,mjx-1))
      allocate(psfc(mix-1,mjx-1))
      allocate(t2(mix-1,mjx-1))
      allocate(q2(mix-1,mjx-1))
      allocate(slp(mix-1,mjx-1))

      allocate(u10response(rix-1,rjx-1))
      allocate(v10response(rix-1,rjx-1))
      allocate(phresponse(rix-1,rjx-1,rkx))
      allocate(phbresponse(rix-1,rjx-1,rkx))
      allocate(tresponse(rix-1,rjx-1,rkx-1))
      allocate(muresponse(rix-1,rjx-1))
      allocate(mubresponse(rix-1,rjx-1))
      allocate(qvaporresponse(rix-1,rjx-1,rkx-1))
      allocate(totpresresponse(rix-1,rjx-1,rkx-1))
      allocate(thetaresponse(rix-1,rjx-1,rkx-1))
      allocate(tempresponse(rix-1,rjx-1,rkx-1))
      allocate(gpttotresponse(rix-1,rjx-1,rkx))
      allocate(gpheightresponse(rix-1,rjx-1,rkx))
      allocate(gpheighthalfresponse(rix-1,rjx-1,rkx-1))
      allocate(rhodryresponse(rix-1,rjx-1,rkx-1))
      allocate(surfptotalresponse(rix-1,rjx-1))
      allocate(psfcresponse(rix-1,rjx-1))
      allocate(t2response(rix-1,rjx-1))
      allocate(q2response(rix-1,rjx-1))
      allocate(slpresponse(rix-1,rjx-1))
      allocate(t850response(rix-1,rjx-1))
      allocate(searainresponsef(rix-1,rjx-1))
      allocate(searainresponsei(rix-1,rjx-1))

      allocate(gph300mean(mix-1,mjx-1))
      allocate(gph500mean(mix-1,mjx-1))
      allocate(gph700mean(mix-1,mjx-1))
      allocate(gph850mean(mix-1,mjx-1))
      allocate(gph925mean(mix-1,mjx-1))
      allocate(u300mean(mix-1,mjx-1))
      allocate(u500mean(mix-1,mjx-1))
      allocate(u700mean(mix-1,mjx-1))
      allocate(u850mean(mix-1,mjx-1))
      allocate(u925mean(mix-1,mjx-1))
      allocate(v300mean(mix-1,mjx-1))
      allocate(v500mean(mix-1,mjx-1))
      allocate(v700mean(mix-1,mjx-1))
      allocate(v850mean(mix-1,mjx-1))
      allocate(v925mean(mix-1,mjx-1))
      allocate(t300mean(mix-1,mjx-1))
      allocate(t500mean(mix-1,mjx-1))
      allocate(t700mean(mix-1,mjx-1))
      allocate(t850mean(mix-1,mjx-1))
      allocate(t925mean(mix-1,mjx-1))
      allocate(gph300(mix-1,mjx-1))
      allocate(gph500(mix-1,mjx-1))
      allocate(gph700(mix-1,mjx-1))
      allocate(gph850(mix-1,mjx-1))
      allocate(gph925(mix-1,mjx-1))
      allocate(u300(mix-1,mjx-1))
      allocate(u500(mix-1,mjx-1))
      allocate(u700(mix-1,mjx-1))
      allocate(u850(mix-1,mjx-1))
      allocate(u925(mix-1,mjx-1))
      allocate(v300(mix-1,mjx-1))
      allocate(v500(mix-1,mjx-1))
      allocate(v700(mix-1,mjx-1))
      allocate(v850(mix-1,mjx-1))
      allocate(v925(mix-1,mjx-1))
      allocate(t300(mix-1,mjx-1))
      allocate(t500(mix-1,mjx-1))
      allocate(t700(mix-1,mjx-1))
      allocate(t850(mix-1,mjx-1))
      allocate(t925(mix-1,mjx-1))
      allocate(uuvar(mix,mjx-1,mkx-1))
      allocate(vvvar(mix-1,mjx,mkx-1))
      allocate(uucrossvar(mix-1,mjx-1,mkx-1))
      allocate(vvcrossvar(mix-1,mjx-1,mkx-1))
      allocate(uucrosscovar(mix-1,mjx-1,mkx-1))
      allocate(vvcrosscovar(mix-1,mjx-1,mkx-1))
      allocate(tempvar(mix-1,mjx-1,mkx-1))
      allocate(totpresvar(mix-1,mjx-1,mkx-1))
      allocate(gph300var(mix-1,mjx-1))
      allocate(gph500var(mix-1,mjx-1))
      allocate(gph700var(mix-1,mjx-1))
      allocate(gph850var(mix-1,mjx-1))
      allocate(gph925var(mix-1,mjx-1))

      allocate(gph300var(mix-1,mjx-1))
      allocate(gph500var(mix-1,mjx-1))
      allocate(gph700var(mix-1,mjx-1))
      allocate(gph850var(mix-1,mjx-1))
      allocate(gph925var(mix-1,mjx-1))

      allocate(uucross(mix-1,mjx-1,mkx-1))
      allocate(vvcross(mix-1,mjx-1,mkx-1))
      allocate(ucrossmean(mix-1,mjx-1,mkx-1))
      allocate(vcrossmean(mix-1,mjx-1,mkx-1))

      allocate(gph300r1covar(mix-1,mjx-1,mkx-1))
      allocate(gph500r1covar(mix-1,mjx-1,mkx-1))
      allocate(gph700r1covar(mix-1,mjx-1,mkx-1))
      allocate(gph850r1covar(mix-1,mjx-1,mkx-1))
      allocate(gph925r1covar(mix-1,mjx-1,mkx-1))
      allocate(gph300targ(mix-1,mjx-1,mkx-1))
      allocate(gph500targ(mix-1,mjx-1,mkx-1))
      allocate(gph700targ(mix-1,mjx-1,mkx-1))
      allocate(gph850targ(mix-1,mjx-1,mkx-1))
      allocate(gph925targ(mix-1,mjx-1,mkx-1))

      allocate(t300r1covar(mix-1,mjx-1,mkx-1))
      allocate(t500r1covar(mix-1,mjx-1,mkx-1))
      allocate(t700r1covar(mix-1,mjx-1,mkx-1))
      allocate(t850r1covar(mix-1,mjx-1,mkx-1))
      allocate(t925r1covar(mix-1,mjx-1,mkx-1))
      allocate(t300targ(mix-1,mjx-1,mkx-1))
      allocate(t500targ(mix-1,mjx-1,mkx-1))
      allocate(t700targ(mix-1,mjx-1,mkx-1))
      allocate(t850targ(mix-1,mjx-1,mkx-1))
      allocate(t925targ(mix-1,mjx-1,mkx-1))

      allocate(u300r1covar(mix-1,mjx-1,mkx-1))
      allocate(u500r1covar(mix-1,mjx-1,mkx-1))
      allocate(u700r1covar(mix-1,mjx-1,mkx-1))
      allocate(u850r1covar(mix-1,mjx-1,mkx-1))
      allocate(u925r1covar(mix-1,mjx-1,mkx-1))
      allocate(u300targ(mix-1,mjx-1,mkx-1))
      allocate(u500targ(mix-1,mjx-1,mkx-1))
      allocate(u700targ(mix-1,mjx-1,mkx-1))
      allocate(u850targ(mix-1,mjx-1,mkx-1))
      allocate(u925targ(mix-1,mjx-1,mkx-1))

      allocate(v300r1covar(mix-1,mjx-1,mkx-1))
      allocate(v500r1covar(mix-1,mjx-1,mkx-1))
      allocate(v700r1covar(mix-1,mjx-1,mkx-1))
      allocate(v850r1covar(mix-1,mjx-1,mkx-1))
      allocate(v925r1covar(mix-1,mjx-1,mkx-1))
      allocate(v300targ(mix-1,mjx-1,mkx-1))
      allocate(v500targ(mix-1,mjx-1,mkx-1))
      allocate(v700targ(mix-1,mjx-1,mkx-1))
      allocate(v850targ(mix-1,mjx-1,mkx-1))
      allocate(v925targ(mix-1,mjx-1,mkx-1))

      allocate(slpr1covar(mix-1,mjx-1,mkx-1))
      allocate(slptarg(mix-1,mjx-1,mkx-1))

      allocate(t2r1covar(mix-1,mjx-1,mkx-1))
      allocate(td2r1covar(mix-1,mjx-1,mkx-1))
      allocate(q2r1covar(mix-1,mjx-1,mkx-1))
      allocate(q850r1covar(mix-1,mjx-1,mkx-1))
      allocate(u10r1covar(mix-1,mjx-1,mkx-1))
      allocate(v10r1covar(mix-1,mjx-1,mkx-1))
      allocate(t2targ(mix-1,mjx-1,mkx-1))
      allocate(td2targ(mix-1,mjx-1,mkx-1))
      allocate(q2targ(mix-1,mjx-1,mkx-1))
      allocate(q850targ(mix-1,mjx-1,mkx-1))
      allocate(u10targ(mix-1,mjx-1,mkx-1))
      allocate(v10targ(mix-1,mjx-1,mkx-1))

      allocate(t300var(mix-1,mjx-1))
      allocate(t500var(mix-1,mjx-1))
      allocate(t700var(mix-1,mjx-1))
      allocate(t850var(mix-1,mjx-1))
      allocate(t925var(mix-1,mjx-1))

      allocate(u300var(mix-1,mjx-1))
      allocate(u500var(mix-1,mjx-1))
      allocate(u700var(mix-1,mjx-1))
      allocate(u850var(mix-1,mjx-1))
      allocate(u925var(mix-1,mjx-1))

      allocate(v300var(mix-1,mjx-1))
      allocate(v500var(mix-1,mjx-1))
      allocate(v700var(mix-1,mjx-1))
      allocate(v850var(mix-1,mjx-1))
      allocate(v925var(mix-1,mjx-1))

      allocate(slpvar(mix-1,mjx-1))
      allocate(psfcvar(mix-1,mjx-1))

      allocate(t2var(mix-1,mjx-1))
      allocate(td2var(mix-1,mjx-1))
      allocate(td2mean(mix-1,mjx-1))
      allocate(td2(mix-1,mjx-1))
      allocate(q2var(mix-1,mjx-1))
      allocate(q850var(mix-1,mjx-1))
      allocate(q850mean(mix-1,mjx-1))
      allocate(q850(mix-1,mjx-1))
      allocate(u10var(mix-1,mjx-1))
      allocate(v10var(mix-1,mjx-1))

      allocate(uggph700(mix-1,mjx-1))
      allocate(uggph850(mix-1,mjx-1))
      allocate(uggph925(mix-1,mjx-1))

      allocate(dbzavgvec(ensnum))
      allocate(dbzmaxvec(ensnum))
      allocate(uhavgvec(ensnum))
      allocate(uhmaxvec(ensnum))
      allocate(pcpvec(ensnum))
      allocate(windavgvec(ensnum))
      allocate(uhcovvec(ensnum))
      allocate(dbzcovvec(ensnum))

c$$$ Response mean allocations

      allocate(dbznext(rix-1,rjx-1,rkx-1))
      allocate(uunext(rix,rjx-1,rkx-1))
      allocate(vvnext(rix-1,rjx,rkx-1))
      allocate(wwnext(rix-1,rjx-1,rkx-1))
      allocate(phnext(rix-1,rjx-1,rkx))
      allocate(phbnext(rix-1,rjx-1,rkx))
      allocate(tnext(rix-1,rjx-1,rkx-1))
      allocate(munext(rix-1,rjx-1))
      allocate(mubnext(rix-1,rjx-1))
      allocate(pnext(rix-1,rjx-1,rkx-1))
      allocate(phydnext(rix-1,rjx-1,rkx-1))
      allocate(pbnext(rix-1,rjx-1,rkx-1))
      allocate(rdnwnext(rkx-1))
      allocate(rdnnext(rkx-1))
      allocate(dnwnext(rkx-1))
      allocate(dnnext(rkx-1))
      allocate(q2next(rix-1,rjx-1))
      allocate(t2next(rix-1,rjx-1))
      allocate(th2next(rix-1,rjx-1))
      allocate(psfcnext(rix-1,rjx-1))
      allocate(u10next(rix-1,rjx-1))
      allocate(v10next(rix-1,rjx-1))
      allocate(qvapornext(rix-1,rjx-1,rkx-1))
      allocate(qrainnext(rix-1,rjx-1,rkx-1))
      allocate(qcloudnext(rix-1,rjx-1,rkx-1))
      allocate(qicenext(rix-1,rjx-1,rkx-1))
      allocate(qsnownext(rix-1,rjx-1,rkx-1))
      allocate(qgraupnext(rix-1,rjx-1,rkx-1))
      allocate(sstnext(rix-1,rjx-1))
      allocate(rainncnext(rix-1,rjx-1))
      allocate(snowncnext(rix-1,rjx-1))
      allocate(hailncnext(rix-1,rjx-1))
      allocate(graupelncnext(rix-1,rjx-1))
      allocate(rainnumnext(rix-1,rjx-1,rkx-1))
      allocate(icenumnext(rix-1,rjx-1,rkx-1))
      allocate(totpresnext(rix-1,rjx-1,rkx-1))
      allocate(uhnext(rix-1,rjx-1))
      allocate(thetanext(rix-1,rjx-1,rkx-1))
      allocate(tempnext(rix-1,rjx-1,rkx-1))
      allocate(gpttotnext(rix-1,rjx-1,rkx))
      allocate(rhodrynext(rix-1,rjx-1,rkx-1))
      allocate(surfptotalnext(rix-1,rjx-1))
      allocate(wspdnext(rix-1,rjx-1))
      allocate(infile(ensnum))

c$$ If specified, save member responses to netcdf file
c$$ First, define response variables

            if (scatter == .true.) then

               filer = 'Rvals.nc'

               rcode = nf_create(filer, nf_clobber, fidr)

               print*, "Created file ", filer

               rcode = nf_def_dim(fidr, 'member', ensnum, mem)

               rcode = nf_def_var(fidr,'DBZ_AVG',nf_float,
     &              1,mem,iddbzavg)
               rcode = nf_put_att_text(fidr,iddbzavg, 'description', 33,
     &              'Average sim. refl. in response function box')
               rcode = nf_put_att_text(fidr,iddbzavg,'units',9,'dBZ')
               rcode = nf_put_att_real(fidr,iddbzavg,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'DBZ_MAX',nf_float,
     &              1,mem,iddbzmax)
               rcode = nf_put_att_text(fidr,iddbzmax, 'description', 33,
     &              'Max sim. refl. in response function box')
               rcode = nf_put_att_text(fidr,iddbzmax,'units',9,'dBZ')
               rcode = nf_put_att_real(fidr,iddbzmax,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'UH_AVG',nf_float,
     &              1,mem,iduhavg)
               rcode = nf_put_att_text(fidr,iduhavg, 'description', 33,
     &              'Average updraft helicity in response function box')
               rcode = nf_put_att_text(fidr,iduhavg,'units',9,
     &              'm**2/s**2')
               rcode = nf_put_att_real(fidr,iduhavg,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'UH_MAX',nf_float,
     &              1,mem,iduhmax)
               rcode = nf_put_att_text(fidr,iduhmax, 'description', 33,
     &              'Max updraft helicity in response function box')
               rcode = nf_put_att_text(fidr,iduhmax,'units',9,
     &              'm**2/s**2')
               rcode = nf_put_att_real(fidr,iduhmax,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'PCP',nf_float,
     &              1,mem,idpcp)
               rcode = nf_put_att_text(fidr,idpcp, 'description', 33,
     &              'Average accum. precip. in response function box')
               rcode = nf_put_att_text(fidr,idpcp,'units',9,
     &              'inches')
               rcode = nf_put_att_real(fidr,idpcp,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'WSPD_AVG',nf_float,
     &              1,mem,idwspd)
               rcode = nf_put_att_text(fidr,idwspd, 'description', 33,
     &              'Average wind speed in response function box')
               rcode = nf_put_att_text(fidr,idwspd,'units',9,
     &              'm/s')
               rcode = nf_put_att_real(fidr,idwspd,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'UH_COV',nf_float,
     &              1,mem,iduhcov)
               rcode = nf_put_att_text(fidr,iduhcov, 'description', 33,
     &              'UH Coverage in response function box')
               rcode = nf_put_att_text(fidr,iduhcov,'units',9,
     &              'm**2/s**2')
               rcode = nf_put_att_real(fidr,iduhcov,'_FillValue',
     &              nf_float,1,fmiss)

               rcode = nf_def_var(fidr,'DBZ_COV',nf_float,
     &              1,mem,iddbzcov)
               rcode = nf_put_att_text(fidr,iddbzcov, 'description', 33,
     &              'Reflectivity Coverage in response function box')
               rcode = nf_put_att_text(fidr,iddbzcov,'units',9,
     &              'dBZ')
               rcode = nf_put_att_real(fidr,iddbzcov,'_FillValue',
     &              nf_float,1,fmiss)

                call close_file(fidr)

            endif

c$$$ Read base state variables from reference file

      call open_file(infilerefvar,permissr,iunitref)
      call get_variable3d(iunitref,'PHB',mix-1,
     &     mjx-1,mkx,1,phbmean)
      call get_variable2d(iunitref,'MUB',mix-1,
     &     mjx-1,1,mub)
      call get_variable2d(iunitref,'MUB',mix-1,
     &     mjx-1,1,mubmean)
      call get_variable3d(iunitref,'PHB',mix-1,
     &     mjx-1,mkx,1,phb)
      call get_variable2d(iunitref,'HGT',mix-1,
     &     mjx-1,1,hgt)
      call close_file(iunitref)

c$$$ Read response means and write to response file
      print*, 'Begin calculating response means'

c$$$ Initialize constants

      outfilemean='Rmean.out'
      iunitwritemean=11

c$$$ Initialize variables for mean calculation
      uhmaxMEM=0.0
      dbzmaxMEM=0.0
      uphelresponse(:,:)=0.0
      dbzresponse(:,:,:)=0.0
      wspresponse(:,:)=0.0
      pcptotresponse(:,:)=0.0

c$$$ Read variables here

c$$$ Beginning of loop through each time
c$$$  to calculate 6-hr response means
c$$$ For speed, only calculating means at
c$$$   response time but can change later
      do cc=1,6

      write(time, '(i2)') cc
      print*, "Time for response mean"
      print*, cc

      print*,infileR(cc, 1)
      call open_file(infileR(cc, 1),permissr,iunitread)
c$$      call get_variable3d(iunitread,'U',rix,
c$$     &     rjx-1,rkx-1,1,uu)
c$$      call get_variable3d(iunitread,'V',rix-1,
c$$     &     rjx,rkx-1,1,vv)
c$$      call get_variable3d(iunitread,'W',rix-1,
c$$     &     rjx-1,rkx-1,1,ww)
c$$      call get_variable3d(iunitread,'PH',rix-1,
c$$     &     rjx-1,rkx,1,ph)
c$$      call get_variable3d(iunitread,'PHB',rix-1,
c$$     &     rjx-1,rkx,0,phb)
c$$      call get_variable3d(iunitread,'T',rix-1,
c$$     &     rjx-1,rkx-1,1,t)
      call get_variable3d(iunitread,'REFL_10CM',rix-1,
     &     rjx-1,rkx-1,1,dbznext)
      call get_variable2d(iunitread,'UP_HELI_MAX',rix-1,
     &     rjx-1,1,uhnext)
c$$      call get_variable2d(iunitread,'MU',rix-1,
c$$     &     rjx-1,1,mu)
c$$      call get_variable2d(iunitread,'MUB',rix-1,
c$$     &     rjx-1,0,mub)
c$$      call get_variable3d(iunitread,'P',rix-1,
c$$     &     rjx-1,rkx-1,1,p)
c$$      call get_variable3d(iunitread,'P_HYD',rix-1,
c$$     &     rjx-1,rkx-1,0,phyd)
c$$      call get_variable3d(iunitread,'PB',rix-1,
c$$     &     rjx-1,rkx-1,0,pb)
c$$      call get_variable2d(iunitread,'Q2',rix-1,
c$$     &     rjx-1,1,q2)
c$$      call get_variable2d(iunitread,'T2',rix-1,
c$$     &     rjx-1,1,t2)
c$$      call get_variable2d(iunitread,'TH2',rix-1,
c$$     &     rjx-1,1,th2)
c$$      call get_variable2d(iunitread,'PSFC',rix-1,
c$$     &     rjx-1,1,psfc)
c$$      call get_variable2d(iunitread,'U10',rix-1,
c$$     &     rjx-1,1,u10)
c$$      call get_variable2d(iunitread,'V10',rix-1,
c$$     &     rjx-1,1,v10)
c$$      call get_variable3d(iunitread,'QVAPOR',rix-1,
c$$     &     rjx-1,rkx-1,1,qvapor)
c$$      call get_variable3d(iunitread,'QRAIN',rix-1,
c$$     &     rjx-1,rkx-1,1,qrain)
c$$      call get_variable3d(iunitread,'QCLOUD',rix-1,
c$$     &     rjx-1,rkx-1,1,qcloud)
c$$      call get_variable3d(iunitread,'QGRAUP',rix-1,
c$$     &     rjx-1,rkx-1,1,qgraup)
c$$      call get_variable3d(iunitread,'QSNOW',rix-1,
c$$     &     rjx-1,rkx-1,1,qsnow)
c$$      call get_variable3d(iunitread,'QICE',rix-1,
c$$     &     rjx-1,rkx-1,1,qice)
c$$      call get_variale2d(iunitread,'SST',rix-1,
c$$     &     rjx-1,0,sst)
      call get_variable2d(iunitread,'RAINNC',rix-1,
     &     rjx-1,1,rainncnext)
      call get_variable2d(iunitread,'WSPD10MAX',rix-1,
     &     rjx-1,1,wspdnext)
c$$      call get_variable2d(iunitread,'SNOWNC',rix-1,
c$$     &     rjx-1,0,snownc)
c$$      call get_variable2d(iunitread,'HAILNC',rix-1,
c$$     &     rjx-1,0,hailnc)
c$$      call get_variable2d(iunitread,'GRAUPELNC',rix-1,
c$$     &     rjx-1,0,graupelnc)
c$$      call get_variable3d(iunitread,'QNRAIN',rix-1,
c$$     &     rjx-1,rkx-1,0,rainnum)
c$$      call get_variable3d(iunitread,'QNICE',rix-1,
c$$     &     rjx-1,rkx-1,0,icenum)
      call close_file(iunitread)

c$$$ Goal is to find the maximum rfunc in rbox for each member over
c$$$  entire six hour time frame.
      if (MAXVAL(uhnext(is:ie,js:je)) .GT. uhmaxMEM) then
          uhmaxMEM = MAXVAL(uhnext(is:ie,js:je))
      end if

      if (MAXVAL(dbznext(is:ie,js:je,1)) .GT. dbzmaxMEM) then
          dbzmaxMEM = MAXVAL(dbznext(is:ie,js:je,1))
      end if

      npts = COUNT((uhnext(is:ie,js:je) .GT. uhthresh))
      uhcovMEM=uhcovMEM+npts

      npts = COUNT((dbznext(is:ie,js:je,1) .GT. dbzthresh))
      dbzcovMEM=dbzcovMEM+npts

      uphelresponse=uphelresponse+uhnext
      dbzresponse=dbzresponse+dbznext
      pcptotresponse=pcptotresponse+rainncnext
      wspresponse=wspresponse+wspdnext

C$$$ End time loop
      enddo

      dbzmaxmean=dbzmaxmean+dbzmaxMEM
      uhmaxmean=uhmaxmean+uhmaxMEM
      uhcovmean=uhcovmean+uhcovMEM
      dbzcovmean=dbzcovmean+dbzcovMEM

      print*,'got all vars for member'
      print*,1


c$$$ Read in all other ens members, summing data for mean
      do m=2,ensnum
      print*, 'Ens mem', m

      uhmaxMEM=0.0
      dbzmaxMEM=0.0
      uhcovMEM=0.0
      dbzcovMEM=0.0

c$$$ Beginning of loop through each time
c$$$  to calculate 6-hr response means
c$$$ For speed, only calculating means at
c$$$   response time but can change later
      do cc=1,6
      write(time, '(i2)') cc
      print*, "Calculating mean data with", infileR(cc, m)

      call open_file(infileR(cc, m),permissr,iunitread)
c$$      call get_variable3d(iunitread,'U',rix,
c$$     &     rjx-1,rkx-1,1,uunext)
c$$      call get_variable3d(iunitread,'V',rix-1,
c$$     &     rjx,rkx-1,1,vvnext)
c$$      call get_variable3d(iunitread,'W',rix-1,
c$$     &     rjx-1,rkx-1,1,wwnext)
c$$      call get_variable3d(iunitread,'PH',rix-1,
c$$     &     rjx-1,rkx,1,phnext)
c$$      call get_variable3d(iunitread,'PHB',rix-1,
c$$     &     rjx-1,rkx,0,phbnext)
c$$      call get_variable3d(iunitread,'T',rix-1,
c$$     &     rjx-1,rkx-1,1,tnext)
      call get_variable3d(iunitread,'REFL_10CM',rix-1,
     &     rjx-1,rkx-1,1,dbznext)
      call get_variable2d(iunitread,'UP_HELI_MAX',rix-1,
     &     rjx-1,1,uhnext)
c$$      call get_variable2d(iunitread,'MU',rix-1,
c$$     &     rjx-1,1,munext)
c$$      call get_variable2d(iunitread,'MUB',rix-1,
c$$     &     rjx-1,0,mubnext)
c$$      call get_variable3d(iunitread,'P',rix-1,
c$$     &    rjx-1,rkx-1,1,pnext)
c$$      call get_variable3d(iunitread,'P_HYD',rix-1,
c$$     &     rjx-1,rkx-1,0,phydnext)
c$$      call get_variable3d(iunitread,'PB',rix-1,
c$$     &     rjx-1,rkx-1,0,pbnext)
c$$      call get_variable1d(iunitread,'RDNW',rkx-1,
c$$     &     1,rdnwnext)
c$$      call get_variable1d(iunitread,'RDN',rkx-1,
c$$     &     1,rdnnext)
c$$      call get_variable2d(iunitread,'Q2',rix-1,
c$$     &     rjx-1,1,q2next)
c$$      call get_variable2d(iunitread,'T2',rix-1,
c$$     &     rjx-1,1,t2next)
c$$      call get_variable2d(iunitread,'TH2',rix-1,
c$$     &     rjx-1,1,th2next)
c$$      call get_variable2d(iunitread,'PSFC',rix-1,
c$$     &     rjx-1,1,psfcnext)
c$$      call get_variable2d(iunitread,'U10',rix-1,
c$$     &     rjx-1,1,u10next)
c$$      call get_variable2d(iunitread,'V10',rix-1,
c$$     &     rjx-1,1,v10next)
c$$      call get_variable3d(iunitread,'QVAPOR',rix-1,
c$$     &     rjx-1,rkx-1,1,qvapornext)
c$$      call get_variable3d(iunitread,'QRAIN',rix-1,
c$$     &     rjx-1,rkx-1,1,qrainnext)
c$$      call get_variable3d(iunitread,'QCLOUD',rix-1,
c$$     &     rjx-1,rkx-1,1,qcloudnext)
c$$      call get_variable3d(iunitread,'QGRAUP',rix-1,
c$$     &     rjx-1,rkx-1,1,qgraupnext)
c$$      call get_variable3d(iunitread,'QSNOW',rix-1,
c$$     &     rjx-1,rkx-1,1,qsnownext)
c$$      call get_variable3d(iunitread,'QICE',rix-1,
c$$     &     rjx-1,rkx-1,1,qicenext)
c$$      call get_variable2d(iunitread,'SST',rix-1,
c$$     &     rjx-1,0,sstnext)
      call get_variable2d(iunitread,'RAINNC',rix-1,
     &     rjx-1,1,rainncnext)
      call get_variable2d(iunitread,'WSPD10MAX',rix-1,
     &     rjx-1,1,wspdnext)
c$$      call get_variable2d(iunitread,'SNOWNC',rix-1,
c$$     &     rjx-1,0,snowncnext)
c$$      call get_variable2d(iunitread,'HAILNC',rix-1,
c$$     &     rjx-1,0,hailncnext)
c$$      call get_variable2d(iunitread,'GRAUPELNC',rix-1,
c$$     &     rjx-1,0,graupelncnext)
c$$      call get_variable3d(iunitread,'QNRAIN',rix-1,
c$$     &     rjx-1,rkx-1,0,rainnumnext)
c$$      call get_variable3d(iunitread,'QNICE',rix-1,
c$$     &     rjx-1,rkx-1,0,icenumnext)

      call close_file(iunitread)

c$$$ Increment max response functions
      if (MAXVAL(uhnext(is:ie,js:je)) .GT. uhmaxMEM) then
          uhmaxMEM = MAXVAL(uhnext(is:ie,js:je))
      end if

      if (MAXVAL(dbznext(is:ie,js:je,1)) .GT. dbzmaxMEM) then
          dbzmaxMEM = MAXVAL(dbznext(is:ie,js:je,1))
      end if

C$$$ Increment coverage response functions
       npts = COUNT((uhnext(is:ie,js:je) .GT. uhthresh))
       print*, 'Number of points in rbox exceeding uh threshold for mem'
       print*, npts
       uhcovMEM=uhcovMEM+npts

       npts = COUNT((dbznext(is:ie,js:je,1) .GT. dbzthresh))
       dbzcovMEM=dbzcovMEM+npts
       print*,'Number of points in rbox exceeding dbz threshold for mem'
       print*, npts

      uphelresponse=uphelresponse+uhnext
      dbzresponse=dbzresponse+dbznext
      pcptotresponse=pcptotresponse+rainncnext
      wspresponse=wspresponse+wspdnext

c$$      uu=uu+uunext
c$$      vv=vv+vvnext
c$$      ww=ww+wwnext
c$$      ph=ph+phnext
c$$      phb=phb+phbnext
c$$      t=t+tnext
c$$      mu=mu+munext
c$$      mub=mub+mubnext
c$$      p=p+pnext
c$$      phyd=phyd+phydnext
c$$      pb=pb+pbnext
c$$      q2=q2+q2next
c$$      t2=t2+t2next
c$$      th2=th2+th2next
c$$      psfc=psfc+psfcnext
c$$      u10=u10+u10next
c$$      v10=v10+v10next
c$$      qvapor=qvapor+qvapornext
c$$      qrain=qrain+qrainnext
c$$      qcloud=qcloud+qcloudnext
c$$      qgraup=qgraup+qgraupnext
c$$      qsnow=qsnow+qsnownext
c$$      qice=qice+qicenext
c$$      sst=sst+sstnext
c$$      snownc=snownc+snowncnext
c$$      rainnum=rainnum+rainnumnext
c$$      icenum=icenum+icenumnext
c$$      hailnc=hailnc+hailncnext
c$$      graupelnc=graupelnc+graupelncnext

c$$$ End of loop through times
      enddo

      dbzmaxmean=dbzmaxmean+dbzmaxMEM
      uhmaxmean=uhmaxmean+uhmaxMEM
      uhcovmean=uhcovmean+uhcovMEM
      dbzcovmean=dbzcovmean+dbzcovMEM

c$$$ End of loop through ens members
      enddo

c$$      uu=uu/ensnum
c$$      vv=vv/ensnum
c$$      ww=ww/ensnum
c$$      ph=ph/ensnum
c$$      phb=phb/ensnum
c$$      t=t/ensnum
c$$      mu=mu/ensnum
c$$      mub=mub/ensnum
c$$      p=p/ensnum
      dbzresponse=dbzresponse/ensnum
      uphelresponse=uphelresponse/ensnum

c$$      phyd=phyd/ensnum
c$$      pb=pb/ensnum
c$$      q2=q2/ensnum
c$$      t2=t2/ensnum
c$$      th2=th2/ensnum
c$$      psfc=psfc/ensnum
c$$      u10=u10/ensnum
c$$      v10=v10/ensnum
c$$      qvapor=qvapor/ensnum
c$$      qrain=qrain/ensnum
c$$      qcloud=qcloud/ensnum
c$$      qsnow=qsnow/ensnum
c$$      qgraup=qgraup/ensnum
c$$      qice=qice/ensnum
c$$      sst=sst/ensnum
      pcptotresponse=pcptotresponse/ensnum
      wspresponse=wspresponse/ensnum
c$$      hailnc=hailnc/ensnum
c$$      snownc=snownc/ensnum
c$$      graupelnc=graupelnc/ensnum
c$$      rainnum=rainnum/ensnum
c$$      icenum=icenum/ensnum

c$$ write response means
      call open_file(outfilemean,permissw,iunitwritemean)

c$$      call write_variable3d(iunitwritemean,'U',rix,
c$$     &     rjx-1,rkx-1,cc,uu)
c$$      call write_variable3d(iunitwritemean,'V',rix-1,
c$$     &     rjx,rkx-1,cc,vv)
c$$      call write_variable3d(iunitwritemean,'W',rix-1,
c$$     &     rjx-1,rkx-1,cc,ww)
c$$      call write_variable3d(iunitwritemean,'PH',rix-1,
c$$     &     rjx-1,rkx,cc,ph)
c$$      call write_variable3d(iunitwritemean,'PHB',mix-1,
c$$     &     mjx-1,mkx,cc,phb)
c$$      call write_variable3d(iunitwritemean,'T',rix-1,
c$$     &     rjx-1,rkx-1,cc,t)
      call write_variable3d(iunitwritemean,'REFL_10CM',rix-1,
     &     rjx-1,rkx-1,timer,dbzresponse)
      call write_variable2d(iunitwritemean,'UP_HELI_MAX',rix-1,
     &     rjx-1,timer,uphelresponse)
c$$      call write_variable2d(iunitwritemean,'MU',rix-1,
c$$     &     rjx-1,cc,mu)
c$$      call write_variable2d(iunitwritemean,'MUB',rix-1,
c$$     &     rjx-1,cc,mub)
c$$      call write_variable3d(iunitwritemean,'P',rix-1,
c$$     &     rjx-1,rkx-1,cc,p)
c$$      call write_variable3d(iunitwritemean,'P_HYD',rix-1,
c$$     &     rjx-1,rkx-1,cc,phyd)
c$$      call write_variable3d(iunitwritemean,'PB',rix-1,
c$$     &     rjx-1,rkx-1,cc,pb)
c$$      call write_variable2d(iunitwritemean,'Q2',rix-1,
c$$     &     rjx-1,cc,q2)
c$$      call write_variable2d(iunitwritemean,'T2',rix-1,
c$$     &     rjx-1,cc,t2)
c$$      call write_variable2d(iunitwritemean,'TH2',rix-1,
c$$     &     rjx-1,cc,th2)
c$$      call write_variable2d(iunitwritemean,'PSFC',rix-1,
c$$     &     rjx-1,cc,psfc)
c$$      call write_variable2d(iunitwritemean,'U10',rix-1,
c$$     &     rjx-1,cc,u10)
c$$      call write_variable2d(iunitwritemean,'V10',rix-1,
c$$     &     rjx-1,cc,v10)
c$$      call write_variable3d(iunitwritemean,'QVAPOR',rix-1,
c$$     &     rjx-1,rkx-1,cc,qvapor)
c$$      call write_variable3d(iunitwritemean,'QRAIN',rix-1,
c$$     &     rjx-1,rkx-1,cc,qrain)
c$$      call write_variable3d(iunitwritemean,'QCLOUD',rix-1,
c$$     &     rjx-1,rkx-1,cc,qcloud)
c$$      call write_variable3d(iunitwritemean,'QSNOW',rix-1,
c$$     &     rjx-1,rkx-1,cc,qsnow)
c$$      call write_variable3d(iunitwritemean,'QGRAUP',rix-1,
c$$     &     rjx-1,rkx-1,cc,qgraup)
c$$      call write_variable3d(iunitwritemean,'QICE',rix-1,
c$$     &     rjx-1,rkx-1,cc,qice)
c$$      call write_variable2d(iunitwritemean,'SST',rix-1,
c$$     &     rjx-1,cc,sst)
      call write_variable2d(iunitwritemean,'RAINNC',rix-1,
     &     rjx-1,timer,pcptotresponse)
      call write_variable2d(iunitwritemean,'WSPD10MAX',rix-1,
     &     rjx-1,timer,wspresponse)
c$$      call write_variable2d(iunitwritemean,'SNOWNC',rix-1,
c$$     &     rjx-1,cc,snownc)
c$$      call write_variable2d(iunitwritemean,'HAILNC',rix-1,
c$$     &     rjx-1,cc,hailnc)
c$$      call write_variable2d(iunitwritemean,'GRAUPELNC',rix-1,
c$$     &     rjx-1,cc,graupelnc)
c$$      call write_variable3d(iunitwritemean,'QNRAIN',rix-1,
c$$     &     rjx-1,rkx-1,cc,rainnum)
c$$      call write_variable3d(iunitwritemean,'QNICE',rix-1,
c$$     &     rjx-1,rkx-1,cc,icenum)

      call close_file(iunitwritemean)

c$$$ Don't divide by six hours because you want mean over
c$$$  entire
      dbzmaxmean=dbzmaxmean/ensnum
      uhmaxmean=uhmaxmean/ensnum
      uhcovmean=uhcovmean/ensnum
      dbzcovmean=dbzcovmean/ensnum

      print*, "DBZ Max Mean: ", dbzmaxmean
      print*, "UH Max Mean: ", uhmaxmean
      print*, "Done with response means"

c$$$ Read sens-time means from file

      print*, "SENS mean file"
      print*, infilemeanSENS
      print*, "SENS time", timesens

      call open_file(infilemeanSENS,permissr,iunitmean)
      call get_variable3d(iunitmean,'U',mix,
     &     mjx-1,mkx-1,timesens,uumean)
      call get_variable3d(iunitmean,'V',mix-1,
     &     mjx,mkx-1,timesens,vvmean)
      call get_variable3d(iunitmean,'PH',mix-1,
     &     mjx-1,mkx,timesens,phmean)
      call get_variable3d(iunitmean,'T',mix-1,
     &     mjx-1,mkx-1,timesens,tmean)
      call get_variable2d(iunitmean,'MU',mix-1,
     &     mjx-1,timesens,mumean)
      call get_variable3d(iunitmean,'QVAPOR',mix-1,
     &     mjx-1,mkx-1,timesens,qvapormean)
      call get_variable2d(iunitmean,'T2',mix-1,
     &     mjx-1,timesens,t2mean)
      call get_variable2d(iunitmean,'RAINNC',mix-1,
     &     mjx-1,timesens,pcptotmean)
c$$      call get_variable2d(iunitmean,'UP_HELI_MAX',mix-1,
c$$     &     mjx-1,timesens,uphelmean)
c$$      call get_variable2d(iunitmean,'WSPD10MAX',mix-1,
c$$     &     mjx-1,timesens,wspmean)
c$$      call get_variable3d(iunitmean,'REFL_10CM',mix-1,
c$$     &     mjx-1,mkx-1,timesens,dbzmean)
c$$      call close_file(iunitmean)

C$$$ Destagger winds

      do k=1,mkx-1
         do i=1,mix-1
            do j=1,mjx-1
               uucrossmean(i,j,k)=(uumean(i,j,k)+
     &              uumean(i+1,j,k))/2
               vvcrossmean(i,j,k)=(vvmean(i,j,k)+
     &              vvmean(i,j+1,k))/2
            enddo
         enddo
      enddo

c$$$ Calculate other variables for sensitivity calcs

c$$$  Transform perturbation theta to temp

      do i=1,mix-1
         do j=1,mjx-1
            do k=1,mkx-1
               thetamean(i,j,k)=tmean(i,j,k)+basethet
            enddo
         enddo
      enddo

c$$$  Calculate total pressure

c$$$  First calculate total geopotential

      do k=1,mkx
         do i=1,mix-1
            do j=1,mjx-1
               gpttotmean(i,j,k)=phmean(i,j,k)+
     &                       phbmean(i,j,k)
            enddo
         enddo
      enddo

c$$$ Calculate geopotential height

      gpheightmean(:,:,:)=gpttotmean(:,:,:)/grav

c$$$  Next calculate dry total surface pressure

      do i=1,mix-1
         do j=1,mjx-1
            surfptotalmean(i,j)=mumean(i,j)+
     &                      mubmean(i,j)
         enddo
      enddo

      do k=1,mkx-1
         levee=znw(k+1)-znw(k)
         do i=1,mix-1
            do j=1,mjx-1
               rhodrymean(i,j,k)=-((levee)/
     &         (gpttotmean(i,j,k+1)-gpttotmean(i,j,k)))*
     &          surfptotalmean(i,j)
            enddo
         enddo
      enddo

c$$$ Now the total pressure
      do i=1,mix-1
         do j=1,mjx-1
            do k=1,mkx-1
               totpresmean(i,j,k)=(thetamean(i,j,k)*prefcont*
     &         rhodrymean(i,j,k)*rdgas*(1.+(rvgas/rdgas)*
     &         qvapormean(i,j,k)))**kapdiv
            enddo
         enddo
      enddo

c$$$ Get geopotential height on half levels

      do i=1,mix-1
         do j=1,mjx-1
            call destag_zstag(znu, znw, mkx-1,
     &           gpheightmean(i,j,:),gpheightmeanhalf(i,j,:))
         enddo
      enddo

c$$$ Interpolate GPH to desired pressure levels

      do i=1,mix-1
         do j=1,mjx-1
            gph300mean(i,j)=interp_pres(gpheightmeanhalf(i,j,:),
     &           totpresmean(i,j,:),30000.,mkx-1)
            gph500mean(i,j)=interp_pres(gpheightmeanhalf(i,j,:),
     &           totpresmean(i,j,:),50000.,mkx-1)
            gph700mean(i,j)=interp_pres(gpheightmeanhalf(i,j,:),
     &           totpresmean(i,j,:),70000.,mkx-1)
            gph850mean(i,j)=interp_pres(gpheightmeanhalf(i,j,:),
     &           totpresmean(i,j,:),85000.,mkx-1)
            gph925mean(i,j)=interp_pres(gpheightmeanhalf(i,j,:),
     &           totpresmean(i,j,:),92500.,mkx-1)
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            u300mean(i,j)=interp_pres(uucrossmean(i,j,:),
     &           totpresmean(i,j,:),30000.,mkx-1)
            u500mean(i,j)=interp_pres(uucrossmean(i,j,:),
     &           totpresmean(i,j,:),50000.,mkx-1)
            u700mean(i,j)=interp_pres(uucrossmean(i,j,:),
     &           totpresmean(i,j,:),70000.,mkx-1)
            u850mean(i,j)=interp_pres(uucrossmean(i,j,:),
     &           totpresmean(i,j,:),85000.,mkx-1)
            u925mean(i,j)=interp_pres(uucrossmean(i,j,:),
     &           totpresmean(i,j,:),92500.,mkx-1)
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            v300mean(i,j)=interp_pres(vvcrossmean(i,j,:),
     &           totpresmean(i,j,:),30000.,mkx-1)
            v500mean(i,j)=interp_pres(vvcrossmean(i,j,:),
     &           totpresmean(i,j,:),50000.,mkx-1)
            v700mean(i,j)=interp_pres(vvcrossmean(i,j,:),
     &           totpresmean(i,j,:),70000.,mkx-1)
            v850mean(i,j)=interp_pres(vvcrossmean(i,j,:),
     &           totpresmean(i,j,:),85000.,mkx-1)
            v925mean(i,j)=interp_pres(vvcrossmean(i,j,:),
     &           totpresmean(i,j,:),92500.,mkx-1)
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            q850mean(i,j)=interp_pres(qvapormean(i,j,:),
     &           totpresmean(i,j,:),85000.,mkx-1)
         enddo
      enddo

c$$$ Now calculate temperature from theta,totpres

      do k=1,mkx-1
         do j=1,mjx-1
            do i=1,mix-1
               tempmean(i,j,k)=thetamean(i,j,k)*
     &         ((preffer/totpresmean(i,j,k))**-kap)
            enddo
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            psfcmean(i,j)=totpresmean(i,j,1)*
     &           (2.71828**((gpheightmeanhalf(i,j,1)
     &           -hgt(i,j))/(29.3*tempmean(i,j,1)
     &           *(1+qvapormean(i,j,1)))))
         enddo
      enddo

c$$$ Calculate SLP

      do i=1,mix-1
         do j=1,mjx-1
            slpmean(i,j)=slp_standard_atmos(t2mean(i,j),
     &           psfcmean(i,j),q2mean(i,j),hgt(i,j))
         enddo
      enddo

c$$$ Calculate dew point

      do i=1,mix-1
         do j=1,mjx-1
            td2mean(i,j)=mixrat_to_tdew(q2mean(i,j),
     &           psfcmean(i,j))
         enddo
      enddo

c$$$ Interpolate U, V, T to desired pressure levels

      do i=1,mix-1
         do j=1,mjx-1
            t300mean(i,j)=interp_pres(tempmean(i,j,:),
     &           totpresmean(i,j,:),30000.,mkx-1)
            t500mean(i,j)=interp_pres(tempmean(i,j,:),
     &           totpresmean(i,j,:),50000.,mkx-1)
            t700mean(i,j)=interp_pres(tempmean(i,j,:),
     &           totpresmean(i,j,:),70000.,mkx-1)
            t850mean(i,j)=interp_pres(tempmean(i,j,:),
     &           totpresmean(i,j,:),85000.,mkx-1)
            t925mean(i,j)=interp_pres(tempmean(i,j,:),
     &           totpresmean(i,j,:),92500.,mkx-1)
         enddo
      enddo

c$$$ Have been having issues with 2m-variable variances
c$$$  so using lowest-sigma level instead

      q2mean = qvapormean(:,:,1)
      u10mean = uucrossmean(:,:,1)
      v10mean = vvcrossmean(:,:,1)

c$$$  Read forecast time means then calculate
c$$$  6-hr response function means

      dbzavg=0.0
      uhavg=0.0
      pcp=0.0
      windavg=0.0
      uhcov=0.0

c$$$ 6-hr means are stored at the response time hour, so only open that
      do cc=timer,timer
          print*, "Re-open response time mean file"
          print*, infilemeanR
          print*, "Current time", cc

          call open_file(infilemeanR,permissr,iunitmean)
          call get_variable3d(iunitmean,'U',rix-1,
     &     rjx-1,1,cc,u10response)
          call get_variable3d(iunitmean,'V',rix-1,
     &     rjx-1,1,cc,v10response)
          call get_variable3d(iunitmean,'PH',rix-1,
     &     rjx-1,rkx,cc,phresponse)
          call get_variable3d(iunitmean,'T',rix-1,
     &     rjx-1,rkx-1,cc,tresponse)
          call get_variable2d(iunitmean,'MU',rix-1,
     &     rjx-1,cc,muresponse)
          call get_variable3d(iunitmean,'QVAPOR',rix-1,
     &     rjx-1,rkx-1,cc,qvaporresponse)
          call get_variable2d(iunitmean,'T2',rix-1,
     &     rjx-1,cc,t2response)
          call get_variable2d(iunitmean,'RAINNC',rix-1,
     &     rjx-1,1,pcptotresponse)
          call get_variable2d(iunitmean,'UP_HELI_MAX',rix-1,
     &     rjx-1,cc,uphelresponse)
          call get_variable2d(iunitmean,'WSPD10MAX',rix-1,
     &     rjx-1,cc,wspresponse)
          call get_variable3d(iunitmean,'REFL_10CM',rix-1,
     &     rjx-1,rkx-1,cc,dbzresponse)
          call close_file(iunitmean)

          q2response = qvaporresponse(:,:,1)

c$$$  Transform perturbation theta to temp

          do i=1,rix-1
             do j=1,rjx-1
                do k=1,rkx-1
                   thetaresponse(i,j,k)=tresponse(i,j,k)+basethet
                enddo
             enddo
          enddo

c$$$  Calculate total pressure

c$$$  First calculate total geopotential

          do k=1,rkx
             do i=1,rix-1
                do j=1,rjx-1
                   gpttotresponse(i,j,k)=phresponse(i,j,k)+
     &                       phbmean(i,j,k)
                enddo
             enddo
          enddo

c$$$  Next calculate dry total surface pressure

          do i=1,rix-1
             do j=1,rjx-1
                surfptotalresponse(i,j)=muresponse(i,j)+
     &                      mubmean(i,j)
             enddo
          enddo

c$$$  Next calculate dry air density

          do k=1,rkx-1
             levee=znw(k+1)-znw(k)
             do i=1,rix-1
                do j=1,rjx-1
                   rhodryresponse(i,j,k)=-((levee)/
     &         (gpttotresponse(i,j,k+1)-gpttotresponse(i,j,k)))*
     &          surfptotalresponse(i,j)
                enddo
             enddo
          enddo

c$$$ Calculate geopotential height

          gpheightresponse(:,:,:)=gpttotresponse(:,:,:)/grav

c$$$ Now the total pressure
          do i=1,rix-1
             do j=1,rjx-1
                do k=1,rkx-1
                   totpresresponse(i,j,k)=(thetaresponse(i,j,k)
     &  *prefcont*rhodryresponse(i,j,k)*rdgas*(1.+(rvgas/rdgas)*
     &         qvaporresponse(i,j,k)))**kapdiv
                enddo
             enddo
          enddo

c$$$ Get geopotential height on half levels

          do i=1,rix-1
             do j=1,rjx-1
                call destag_zstag(znu, znw, mkx-1,
     &           gpheightresponse(i,j,:),
     &           gpheighthalfresponse(i,j,:))
             enddo
          enddo

c$$$ Now calculate temperature from theta,totpres

          do k=1,rkx-1
             do j=1,rjx-1
                do i=1,rix-1
                   tempresponse(i,j,k)=thetaresponse(i,j,k)*
     &         ((preffer/totpresresponse(i,j,k))**-kap)
                enddo
             enddo
          enddo

c$$$ Calculate Surface Pressure

          do i=1,rix-1
             do j=1,rjx-1
                psfcresponse(i,j)=totpresresponse(i,j,1)*
     &           (2.71828**((gpheighthalfresponse(i,j,1)
     &           -hgt(i,j))/(29.3*tempresponse(i,j,1)
     &           *(1+qvaporresponse(i,j,1)))))
             enddo
          enddo

c$$$ Calculate SLP

          do i=1,rix-1
             do j=1,rjx-1
                slpresponse(i,j)=slp_standard_atmos(t2response(i,j),
     &           psfcresponse(i,j),q2response(i,j),hgt(i,j))
             enddo
          enddo

C$$$ Calculate Response function means of choice

C$$$ Increment 6-h Avg sim dBZ

          do i=is,ie
             do j=js,je
                         dbzavg=dbzavg+dbzresponse(i,j,1)
             enddo
          enddo

C$$$ Inrement uphel

          do i=is,ie
             do j=js,je
                uhavg=uhavg+uphelresponse(i,j)
             enddo
          enddo

C$$$ Accum pcp

          do i=is,ie
             do j=js,je
                pcp=pcp+pcptotresponse(i,j)
             enddo
          enddo


C$$$ Avg max 10-m winds

          do i=is,ie
             do j=js,je
                windavg=windavg+wspresponse(i,j)
             enddo
          enddo
c$$ end of 6-hr time frame
      enddo

c$$ Finish response box means
      dbzavg=dbzavg/((ie-is+1)*(je-js+1)*6)
      uhavg=uhavg/((ie-is+1)*(je-js+1)*6)
      pcp=(pcp/((ie-is+1)*(je-js+1)*6))/25.4
      windavg=windavg/((ie-is+1)*(je-js+1)*6)

C$$$ Need variance of all sens variables

      gph300var(:,:)=0.0
      gph500var(:,:)=0.0
      gph700var(:,:)=0.0
      gph850var(:,:)=0.0
      gph925var(:,:)=0.0
      t300var(:,:)=0.0
      t500var(:,:)=0.0
      t700var(:,:)=0.0
      t850var(:,:)=0.0
      t925var(:,:)=0.0
      v300var(:,:)=0.0
      v500var(:,:)=0.0
      v700var(:,:)=0.0
      v850var(:,:)=0.0
      v925var(:,:)=0.0
      u300var(:,:)=0.0
      u500var(:,:)=0.0
      u700var(:,:)=0.0
      u850var(:,:)=0.0
      u925var(:,:)=0.0
      slpvar(:,:)=0.0
      t2var(:,:)=0.0
      td2var(:,:)=0.0
      q2var(:,:)=0.0
      q850var(:,:)=0.0
      u10var(:,:)=0.0
      v10var(:,:)=0.0

C$$$ Need covar

      gph300r1covar(:,:,:)=0.0
      gph500r1covar(:,:,:)=0.0
      gph700r1covar(:,:,:)=0.0
      gph850r1covar(:,:,:)=0.0
      gph925r1covar(:,:,:)=0.0

      v300r1covar(:,:,:)=0.0
      v500r1covar(:,:,:)=0.0
      v700r1covar(:,:,:)=0.0
      v850r1covar(:,:,:)=0.0
      v925r1covar(:,:,:)=0.0

      u300r1covar(:,:,:)=0.0
      u500r1covar(:,:,:)=0.0
      u700r1covar(:,:,:)=0.0
      u850r1covar(:,:,:)=0.0
      u925r1covar(:,:,:)=0.0

      t300r1covar(:,:,:)=0.0
      t500r1covar(:,:,:)=0.0
      t700r1covar(:,:,:)=0.0
      t850r1covar(:,:,:)=0.0
      t925r1covar(:,:,:)=0.0

      slpr1covar(:,:,:)=0.0
      t2r1covar(:,:,:)=0.0
      td2r1covar(:,:,:)=0.0
      q2r1covar(:,:,:)=0.0
      q850r1covar(:,:,:)=0.0
      u10r1covar(:,:,:)=0.0
      v10r1covar(:,:,:)=0.0

      uggph700(:,:)=0
      uggph850(:,:)=0
      uggph925(:,:)=0

      smat(:,:,:)=0.0
      tmat(:,:,:)=0.0
      mmat(:,:,:)=0.0

c$$$ Read in all ens members

      do m=1,ensnum

c$$$ first sens-time vars

      print*, 'Current SENS Member: ' // infileSENS(m)

      call open_file(infileSENS(m),permissr,iunitread)
      call get_variable3d(iunitread,'U',mix,
     &     mjx-1,mkx-1,1,uu)
      call get_variable3d(iunitread,'V',mix-1,
     &     mjx,mkx-1,1,vv)
      call get_variable3d(iunitread,'PH',mix-1,
     &     mjx-1,mkx,1,ph)
      call get_variable3d(iunitread,'T',mix-1,
     &     mjx-1,mkx-1,1,t)
      call get_variable2d(iunitread,'MU',mix-1,
     &     mjx-1,1,mu)
      call get_variable3d(iunitread,'QVAPOR',mix-1,
     &     mjx-1,mkx-1,1,qvapor)
      call get_variable2d(iunitread,'T2',mix-1,
     &     mjx-1,1,t2)
      call get_variable2d(iunitread,'RAINNC',mix-1,
     &     mjx-1,1,pcptot)
      call get_variable2d(iunitread,'UP_HELI_MAX',mix-1,
     &     mjx-1,1,uphel)
      call get_variable2d(iunitread,'WSPD10MAX',mix-1,
     &     mjx-1,1,wsp)
      call get_variable3d(iunitread,'REFL_10CM',mix-1,
     &     mjx-1,mkx-1,1,dbz)
      call close_file(iunitread)

C$$$ Destagger winds

      do k=1,mkx-1
         do i=1,mix-1
            do j=1,mjx-1
               uucross(i,j,k)=(uu(i,j,k)+
     &              uu(i+1,j,k))/2
               vvcross(i,j,k)=(vv(i,j,k)+
     &              vv(i,j+1,k))/2
            enddo
         enddo
      enddo

c$$$ Calculate other variables for sensitivity calcs

c$$$  Transform perturbation theta to temp

      do i=1,mix-1
         do j=1,mjx-1
            do k=1,mkx-1
               theta(i,j,k)=t(i,j,k)+basethet
            enddo
         enddo
      enddo

c$$$  Calculate total pressure

c$$$  First calculate total geopotential

      do k=1,mkx
         do i=1,mix-1
            do j=1,mjx-1
               gpttot(i,j,k)=ph(i,j,k)+
     &                       phb(i,j,k)
            enddo
         enddo
      enddo

c$$$ Calculate geopotential height

      gpheight(:,:,:)=gpttot(:,:,:)/grav

c$$$  Next calculate dry total surface pressure

      do i=1,mix-1
         do j=1,mjx-1
            surfptotal(i,j)=mu(i,j)+
     &                      mub(i,j)
         enddo
      enddo

      do k=1,mkx-1
         levee=znw(k+1)-znw(k)
         do i=1,mix-1
            do j=1,mjx-1
               rhodry(i,j,k)=-((levee)/
     &         (gpttot(i,j,k+1)-gpttot(i,j,k)))*
     &          surfptotal(i,j)
            enddo
         enddo
      enddo

c$$$ Now the total pressure
      do i=1,mix-1
         do j=1,mjx-1
            do k=1,mkx-1
               totpres(i,j,k)=(theta(i,j,k)*prefcont*
     &         rhodry(i,j,k)*rdgas*(1.+(rvgas/rdgas)*
     &         qvapor(i,j,k)))**kapdiv
            enddo
         enddo
      enddo

c$$$ Get geopotential height on half levels

      do i=1,mix-1
         do j=1,mjx-1
            call destag_zstag(znu, znw, mkx-1,
     &           gpheight(i,j,:),gpheighthalf(i,j,:))
         enddo
      enddo

c$$$ Interpolate GPH to desired pressure levels

      do i=1,mix-1
         do j=1,mjx-1
            gph300(i,j)=interp_pres(gpheighthalf(i,j,:),
     &           totpres(i,j,:),30000.,mkx-1)
            gph500(i,j)=interp_pres(gpheighthalf(i,j,:),
     &           totpres(i,j,:),50000.,mkx-1)
            gph700(i,j)=interp_pres(gpheighthalf(i,j,:),
     &           totpres(i,j,:),70000.,mkx-1)
            gph850(i,j)=interp_pres(gpheighthalf(i,j,:),
     &           totpres(i,j,:),85000.,mkx-1)
            gph925(i,j)=interp_pres(gpheighthalf(i,j,:),
     &           totpres(i,j,:),92500.,mkx-1)
         enddo
      enddo

c$$$ Now calculate temperature from theta,totpres

      do k=1,mkx-1
         do j=1,mjx-1
            do i=1,mix-1
               temp(i,j,k)=theta(i,j,k)*
     &         ((preffer/totpres(i,j,k))**-kap)
            enddo
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            psfc(i,j)=totpres(i,j,1)*
     &           (2.71828**((gpheighthalf(i,j,1)
     &           -hgt(i,j))/(29.3*temp(i,j,1)
     &           *(1+qvapor(i,j,1)))))
         enddo
      enddo

c$$$ Calculate SLP

      do i=1,mix-1
         do j=1,mjx-1
            slp(i,j)=slp_standard_atmos(t2(i,j),
     &           psfc(i,j),q2(i,j),hgt(i,j))
         enddo
      enddo

c$$$ Calculate dew point

      do i=1,mix-1
         do j=1,mjx-1
            td2(i,j)=mixrat_to_tdew(q2(i,j),
     &           psfc(i,j))
         enddo
      enddo

c$$$ Interpolate U, V, T to desired pressure levels

      do i=1,mix-1
         do j=1,mjx-1
            t300(i,j)=interp_pres(temp(i,j,:),
     &           totpres(i,j,:),30000.,mkx-1)
            t500(i,j)=interp_pres(temp(i,j,:),
     &           totpres(i,j,:),50000.,mkx-1)
            t700(i,j)=interp_pres(temp(i,j,:),
     &           totpres(i,j,:),70000.,mkx-1)
            t850(i,j)=interp_pres(temp(i,j,:),
     &           totpres(i,j,:),85000.,mkx-1)
            t925(i,j)=interp_pres(temp(i,j,:),
     &           totpres(i,j,:),92500.,mkx-1)
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            u300(i,j)=interp_pres(uucross(i,j,:),
     &           totpres(i,j,:),30000.,mkx-1)
            u500(i,j)=interp_pres(uucross(i,j,:),
     &           totpres(i,j,:),50000.,mkx-1)
            u700(i,j)=interp_pres(uucross(i,j,:),
     &           totpres(i,j,:),70000.,mkx-1)
            u850(i,j)=interp_pres(uucross(i,j,:),
     &           totpres(i,j,:),85000.,mkx-1)
            u925(i,j)=interp_pres(uucross(i,j,:),
     &           totpres(i,j,:),92500.,mkx-1)
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            v300(i,j)=interp_pres(vvcross(i,j,:),
     &           totpres(i,j,:),30000.,mkx-1)
            v500(i,j)=interp_pres(vvcross(i,j,:),
     &           totpres(i,j,:),50000.,mkx-1)
            v700(i,j)=interp_pres(vvcross(i,j,:),
     &           totpres(i,j,:),70000.,mkx-1)
            v850(i,j)=interp_pres(vvcross(i,j,:),
     &           totpres(i,j,:),85000.,mkx-1)
            v925(i,j)=interp_pres(vvcross(i,j,:),
     &           totpres(i,j,:),92500.,mkx-1)
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            q850(i,j)=interp_pres(qvapor(i,j,:),
     &           totpres(i,j,:),85000.,mkx-1)
         enddo
      enddo

c$$$ Check for gph underground
      do i=1,mix-1
         do j=1,mjx-1
            if (gph700(i,j) .LT. 0.0) then
               uggph700(i,j)=uggph700(i,j)+1
            endif
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            if (gph850(i,j) .LT. 0.0) then
               uggph850(i,j)=uggph850(i,j)+1
            endif
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            if (gph925(i,j) .LT. 0.0) then
               uggph925(i,j)=uggph925(i,j)+1
            endif
         enddo
      enddo

c$$$      do i=1,mix-1
c$$$         do j=1,mjx-1
c$$$            do k=1,mkx-1
c$$$               uucrossvar(i,j,k)=uucrossvar(i,j,k)+
c$$$     &         (uucross(i,j,k)-uucrossmean(i,j,k))**2
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$      do i=1,mix-1
c$$$         do j=1,mjx-1
c$$$            do k=1,mkx-1
c$$$               vvcrossvar(i,j,k)=vvcrossvar(i,j,k)+
c$$$     &         (vvcross(i,j,k)-vvcrossmean(i,j,k))**2
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$      do i=1,mix-1
c$$$         do j=1,mjx-1
c$$$            do k=1,mkx-1
c$$$               tempvar(i,j,k)=tempvar(i,j,k)+
c$$$     &         (temp(i,j,k)-tempmean(i,j,k))**2
c$$$            enddo
c$$$         enddo
c$$$      enddo

c$$$      do i=1,mix-1
c$$$         do j=1,mjx-1
c$$$            do k=1,mkx-1
c$$$               totpresvar(i,j,k)=totpresvar(i,j,k)+
c$$$     &         (totpres(i,j,k)-totpresmean(i,j,k))**2
c$$$            enddo
c$$$         enddo
c$$$      enddo

      do i=1,mix-1
         do j=1,mjx-1
            slpvar(i,j)=slpvar(i,j)+
     &      (slp(i,j)-slpmean(i,j))**2
         enddo
      enddo

c$$$      do i=1,mix-1
c$$$         do j=1,mjx-1
c$$$            psfcvar(i,j)=psfcvar(i,j)+
c$$$     &      (psfc(i,j)-psfcmean(i,j))**2
c$$$         enddo
c$$$      enddo


      do i=1,mix-1
         do j=1,mjx-1
            gph300var(i,j)=gph300var(i,j)+
     &      (gph300(i,j)-gph300mean(i,j))**2
            gph500var(i,j)=gph500var(i,j)+
     &      (gph500(i,j)-gph500mean(i,j))**2
            gph700var(i,j)=gph700var(i,j)+
     &      (gph700(i,j)-gph700mean(i,j))**2
            gph850var(i,j)=gph850var(i,j)+
     &      (gph850(i,j)-gph850mean(i,j))**2
            gph925var(i,j)=gph925var(i,j)+
     &      (gph925(i,j)-gph925mean(i,j))**2
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            t300var(i,j)=t300var(i,j)+
     &      (t300(i,j)-t300mean(i,j))**2
            t500var(i,j)=t500var(i,j)+
     &      (t500(i,j)-t500mean(i,j))**2
            t700var(i,j)=t700var(i,j)+
     &      (t700(i,j)-t700mean(i,j))**2
            t850var(i,j)=t850var(i,j)+
     &      (t850(i,j)-t850mean(i,j))**2
            t925var(i,j)=t925var(i,j)+
     &      (t925(i,j)-t925mean(i,j))**2
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            u300var(i,j)=u300var(i,j)+
     &      (u300(i,j)-u300mean(i,j))**2
            u500var(i,j)=u500var(i,j)+
     &      (u500(i,j)-u500mean(i,j))**2
            u700var(i,j)=u700var(i,j)+
     &      (u700(i,j)-u700mean(i,j))**2
            u850var(i,j)=u850var(i,j)+
     &      (u850(i,j)-u850mean(i,j))**2
            u925var(i,j)=u925var(i,j)+
     &      (u925(i,j)-u925mean(i,j))**2
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            v300var(i,j)=v300var(i,j)+
     &      (v300(i,j)-v300mean(i,j))**2
            v500var(i,j)=v500var(i,j)+
     &      (v500(i,j)-v500mean(i,j))**2
            v700var(i,j)=v700var(i,j)+
     &      (v700(i,j)-v700mean(i,j))**2
            v850var(i,j)=v850var(i,j)+
     &      (v850(i,j)-v850mean(i,j))**2
            v925var(i,j)=v925var(i,j)+
     &      (v925(i,j)-v925mean(i,j))**2
         enddo
      enddo

      do i=1,mix-1
         do j=1,mjx-1
            q850var(i,j)=q850var(i,j)+
     &      (q850(i,j)-q850mean(i,j))**2
            t2var(i,j)=t2var(i,j)+
     &      (t2(i,j)-t2mean(i,j))**2
            q2var(i,j)=q2var(i,j)+
     &      (q2(i,j)-q2mean(i,j))**2
            td2var(i,j)=td2var(i,j)+
     &      (td2(i,j)-td2mean(i,j))**2
            u10var(i,j)=u10var(i,j)+
     &      (u10(i,j)-u10mean(i,j))**2
            v10var(i,j)=v10var(i,j)+
     &      (v10(i,j)-v10mean(i,j))**2
         enddo
      enddo
      q2 = qvapor(:,:,1)
      u10 = uucross(:,:,1)
      v10 = vvcross(:,:,1)

c      print*, "U10 Lowest variance", u10(u10maxii, u10maxjj)
      print*, "V10 variance", v10var(v10maxii, v10maxjj)
      print*, "V10 Mean", v10mean(v10maxii, v10maxjj)
      print*, "V10 Lowest variance", v10(v10maxii, v10maxjj)

c      print*, "Q2 Lowest variance", q2(q2maxii, q2maxjj)

c$$$ Now get all forecast-time response funciton vars

c$$$ Initialize member values
      dbzavgMEM=0.0
      uhavgMEM=0.0
      pcpMEM=0.0
      windavgMEM=0.0
      uhmaxMEM=0.0
      dbzmaxMEM=0.0
      uhcovMEM=0.0
      dbzcovMEM=0.0

      do cc=1,6
      print*, 'Current Response Member: ' // infileR(cc,m)

      call open_file(infileR(cc,m),permissr,iunitread)
      call get_variable3d(iunitread,'PH',rix-1,
     &     rjx-1,rkx,1,phresponse)
      call get_variable3d(iunitread,'T',rix-1,
     &     rjx-1,rkx-1,1,tresponse)
      call get_variable2d(iunitread,'MU',rix-1,
     &     rjx-1,1,muresponse)
      call get_variable3d(iunitread,'QVAPOR',rix-1,
     &     rjx-1,rkx-1,1,qvaporresponse)
      call get_variable2d(iunitread,'T2',rix-1,
     &     rjx-1,1,t2response)
      call get_variable2d(iunitread,'RAINNC',rix-1,
     &     rjx-1,1,pcptotresponse)
      call get_variable2d(iunitread,'UP_HELI_MAX',rix-1,
     &     rjx-1,1,uphelresponse)
      call get_variable2d(iunitread,'WSPD10MAX',rix-1,
     &     rjx-1,1,wspresponse)
      call get_variable3d(iunitread,'REFL_10CM',rix-1,
     &     rjx-1,rkx-1,1,dbzresponse)
      call close_file(iunitread)
      q2response = qvaporresponse(:,:,1)

c$$$  Transform perturbation theta to temp

      do i=1,rix-1
         do j=1,rjx-1
            do k=1,rkx-1
               thetaresponse(i,j,k)=tresponse(i,j,k)+basethet
            enddo
         enddo
      enddo

c$$$  Calculate total pressure

c$$$  First calculate total geopotential

      do k=1,rkx
         do i=1,rix-1
            do j=1,rjx-1
               gpttotresponse(i,j,k)=phresponse(i,j,k)+
     &                       phbmean(i,j,k)
            enddo
         enddo
      enddo

c$$$  Next calculate dry total surface pressure

      do i=1,rix-1
         do j=1,rjx-1
            surfptotalresponse(i,j)=muresponse(i,j)+
     &                      mubmean(i,j)
         enddo
      enddo

c$$$  Next calculate dry air density

      do k=1,rkx-1
         levee=znw(k+1)-znw(k)
         do i=1,rix-1
            do j=1,rjx-1
               rhodryresponse(i,j,k)=-((levee)/
     &         (gpttotresponse(i,j,k+1)-gpttotresponse(i,j,k)))*
     &          surfptotalresponse(i,j)
            enddo
         enddo
      enddo

c$$$ Calculate geopotential height

      gpheightresponse(:,:,:)=gpttotresponse(:,:,:)/grav

c$$$ Now the total pressure
      do i=1,rix-1
         do j=1,rjx-1
            do k=1,rkx-1
               totpresresponse(i,j,k)=(thetaresponse(i,j,k)
     &  *prefcont*rhodryresponse(i,j,k)*rdgas*(1.+(rvgas/rdgas)*
     &         qvaporresponse(i,j,k)))**kapdiv
            enddo
         enddo
      enddo

c$$$ Get geopotential height on half levels

      do i=1,rix-1
         do j=1,rjx-1
            call destag_zstag(znu, znw, mkx-1,
     &           gpheightresponse(i,j,:),
     &           gpheighthalfresponse(i,j,:))
         enddo
      enddo

c$$$ Now calculate temperature from theta,totpres

      do k=1,rkx-1
         do j=1,rjx-1
            do i=1,rix-1
               tempresponse(i,j,k)=thetaresponse(i,j,k)*
     &         ((preffer/totpresresponse(i,j,k))**-kap)
            enddo
         enddo
      enddo

c$$$ Calculate Surface Pressure

      do i=1,rix-1
         do j=1,rjx-1
            psfcresponse(i,j)=totpresresponse(i,j,1)*
     &           (2.71828**((gpheighthalfresponse(i,j,1)
     &           -hgt(i,j))/(29.3*tempresponse(i,j,1)
     &           *(1+qvaporresponse(i,j,1)))))
         enddo
      enddo

c$$$ Calculate SLP

      do i=1,rix-1
         do j=1,rjx-1
            slpresponse(i,j)=slp_standard_atmos(t2response(i,j),
     &           psfcresponse(i,j),q2response(i,j),hgt(i,j))
         enddo
      enddo

C$$$ Calculate R of choice

C$$$ Avg sim dBZ


      do i=is,ie
         do j=js,je
            dbzavgMEM=dbzavgMEM+dbzresponse(i,j,1)
         enddo
      enddo

C$$$ Max dBZ

      if (MAXVAL(dbzresponse(is:ie,js:je,1)) > dbzmaxMEM) then
          dbzmaxMEM=MAXVAL(dbzresponse(is:ie,js:je,1))
      end if

C$$$ Avg uphel

      do i=is,ie
          do j=js,je
              uhavgMEM=uhavgMEM+uphelresponse(i,j)
          enddo
      enddo

C$$$ Max uphel

      if (MAXVAL(uphelresponse(is:ie,js:je)) > uhmaxMEM) then
          uhmaxMEM=MAXVAL(uphelresponse(is:ie,js:je))
      end if

C$$$ Accum pcp

      do i=is,ie
         do j=js,je
            pcpMEM=pcpMEM+pcptotresponse(i,j)
         enddo
      enddo

C$$$ Avg 10-m winds

      do i=is,ie
         do j=js,je
            windavgMEM=windavgMEM+wspresponse(i,j)
         enddo
      enddo

C$$$ Increment coverage response functions

c$$$ UH Coverage
      npts = COUNT((uphelresponse(is:ie,js:je) .GT. uhthresh))
      print*, "Number of points in rbox exceeding uh threshold for mem"
      print*, npts
      uhcovMEM=uhcovMEM+npts

c$$$ DBZ Coverage
      npts = COUNT((dbzresponse(is:ie,js:je,1) .GT. dbzthresh))
      dbzcovMEM=dbzcovMEM+npts
      print*,"Number of points in rbox exceeding dbz threshold for mem"
      print*, npts

c$$$ End 6-hr time loop
      enddo

c$$$ Finish six-hour member averages
      dbzavgMEM=dbzavgMEM/((ie-is+1)*(je-js+1)*6)
      uhavgMEM=uhavgMEM/((ie-is+1)*(je-js+1)*6)
      pcpMEM=(pcpMEM/((ie-is+1)*(je-js+1))*6)/25.4
      windavgMEM=windavgMEM/((ie-is+1)*(je-js+1)*6)

c$$$ Assemble R perts from mean in column vector

      print*, "Avg dBZ for member", dbzavgMEM
      print*, "Mean Avg dBZ", dbzavg
      print*, "Max dBZ for member", dbzmaxMEM
      print*, "Mean Max dBZ", dbzmaxmean
      print*, "Avg UH for member", uhavgMEM
      print*, "Mean Avg UH", uhavg
      print*, "Max UH for member", uhmaxMEM
      print*, "Mean Max UH", uhmaxmean
      print*, "Avg Accum Precip for member", pcpMEM
      print*, "Mean Avg Accump Precip", pcp
      print*, "Avg Wind for member", windavgMEM
      print*, "Mean Avg Wind", windavg
      print*, "UH Coverage for member", uhcovMEM
      print*, "Mean UH Coverage", uhcovmean
      print*, "DBZ Coverage for member", dbzcovMEM
      print*, "Mean DBZ Coverage", dbzcovmean

      evec1(1)=dbzavgMEM-dbzavg
      evec1(2)=dbzmaxMEM-dbzmaxmean
      evec1(3)=uhavgMEM-uhavg
      evec1(4)=uhmaxMEM-uhmaxmean
      evec1(5)=pcpMEM-pcp
      evec1(6)=windavgMEM-windavg
      evec1(7)=uhcovMEM-uhcovmean
      evec1(8)=dbzcovMEM-dbzcovmean

c$$$ Create vector of responses

      dbzavgvec(m) = dbzavgMEM
      dbzmaxvec(m) = dbzmaxMEM
      uhavgvec(m) = uhavgMEM
      uhmaxvec(m) = uhmaxMEM
      pcpvec(m) = pcpMEM
      windavgvec(m) = windavgMEM
      uhcovvec(m) = uhcovMEM
      dbzcovvec(m) = dbzcovMEM

c$$$ Calculate covariances b/w response function, IC vars
      print*, evec1(rrr)
      do i=1,mix-1
         do j=1,mjx-1
            do RR=rrr,rrr
            gph300r1covar(i,j,RR)=gph300r1covar(i,j,RR)+
     &      ((gph300(i,j)-gph300mean(i,j))*evec1(RR))
            gph500r1covar(i,j,RR)=gph500r1covar(i,j,RR)+
     &      ((gph500(i,j)-gph500mean(i,j))*evec1(RR))
            gph700r1covar(i,j,RR)=gph700r1covar(i,j,RR)+
     &      ((gph700(i,j)-gph700mean(i,j)))*evec1(RR)
            gph850r1covar(i,j,RR)=gph850r1covar(i,j,RR)+
     &      ((gph850(i,j)-gph850mean(i,j)))*evec1(RR)
            gph925r1covar(i,j,RR)=gph925r1covar(i,j,RR)+
     &      ((gph925(i,j)-gph925mean(i,j)))*evec1(RR)


            t300r1covar(i,j,RR)=t300r1covar(i,j,RR)+
     &      (t300(i,j)-t300mean(i,j))*evec1(RR)
            t500r1covar(i,j,RR)=t500r1covar(i,j,RR)+
     &      (t500(i,j)-t500mean(i,j))*evec1(RR)
            t700r1covar(i,j,RR)=t700r1covar(i,j,RR)+
     &      (t700(i,j)-t700mean(i,j))*evec1(RR)
            t850r1covar(i,j,RR)=t850r1covar(i,j,RR)+
     &      (t850(i,j)-t850mean(i,j))*evec1(RR)
            t925r1covar(i,j,RR)=t925r1covar(i,j,RR)+
     &      (t925(i,j)-t925mean(i,j))*evec1(RR)


            u300r1covar(i,j,RR)=u300r1covar(i,j,RR)+
     &      (u300(i,j)-u300mean(i,j))*evec1(RR)
            u500r1covar(i,j,RR)=u500r1covar(i,j,RR)+
     &      (u500(i,j)-u500mean(i,j))*evec1(RR)
            u700r1covar(i,j,RR)=u700r1covar(i,j,RR)+
     &      (u700(i,j)-u700mean(i,j))*evec1(RR)
            u850r1covar(i,j,RR)=u850r1covar(i,j,RR)+
     &      (u850(i,j)-u850mean(i,j))*evec1(RR)
            u925r1covar(i,j,RR)=u925r1covar(i,j,RR)+
     &      (u925(i,j)-u925mean(i,j))*evec1(RR)

            v300r1covar(i,j,RR)=v300r1covar(i,j,RR)+
     &      (v300(i,j)-v300mean(i,j))*evec1(RR)
            v500r1covar(i,j,RR)=v500r1covar(i,j,RR)+
     &      (v500(i,j)-v500mean(i,j))*evec1(RR)
            v700r1covar(i,j,RR)=v700r1covar(i,j,RR)+
     &      (v700(i,j)-v700mean(i,j))*evec1(RR)
            v850r1covar(i,j,RR)=v850r1covar(i,j,RR)+
     &      (v850(i,j)-v850mean(i,j))*evec1(RR)
            v925r1covar(i,j,RR)=v925r1covar(i,j,RR)+
     &      (v925(i,j)-v925mean(i,j))*evec1(RR)

            slpr1covar(i,j,RR)=slpr1covar(i,j,RR)+
     &      (slp(i,j)-slpmean(i,j))*evec1(RR)

            q850r1covar(i,j,RR)=q850r1covar(i,j,RR)+
     &      (q850(i,j)-q850mean(i,j))*evec1(RR)

            t2r1covar(i,j,RR)=t2r1covar(i,j,RR)+
     &      (t2(i,j)-t2mean(i,j))*evec1(RR)

            q2r1covar(i,j,RR)=q2r1covar(i,j,RR)+
     &      (q2(i,j)-q2mean(i,j))*evec1(RR)

            td2r1covar(i,j,RR)=td2r1covar(i,j,RR)+
     &      (td2(i,j)-td2mean(i,j))*evec1(RR)

            v10r1covar(i,j,RR)=v10r1covar(i,j,RR)+
     &      (v10(i,j)-v10mean(i,j))*evec1(RR)

            u10r1covar(i,j,RR)=u10r1covar(i,j,RR)+
     &      (u10(i,j)-u10mean(i,j))*evec1(RR)

            enddo

         enddo
      enddo

      print*,'got all vars for member'
      print*,m

c$$$ End of loop through ens members
      enddo

        do i=1,mix-1
        do j=1,mjx-1
           if (q2var(i,j) .eq. minval(q2var)) then
              print*, "MIN Q2 VARIANCE", minval(q2var)
              q2maxii = i
              q2maxjj = j
           endif
           if (u10var(i,j) .eq. minval(u10var)) then
              print*, "MIN U10 VARIANCE", minval(u10var)
              u10maxii = i
              u10maxjj = j
           endif
           if (v10var(i,j) .eq. minval(v10var)) then
              print*, "MIN V10 VARIANCE", minval(v10var)
              v10maxii = i
              v10maxjj = j
           endif
        enddo
       enddo
      print*, "Q2 Min var loc", q2maxii, q2maxjj
      print*, "V10 min var loc", v10maxii, v10maxjj
      print*, "U10 min var loc", u10maxii, u10maxjj

C$$$ Calcualte esens, store in covar variables

c$$$      do RR=1,numresp
      do RR=rrr,rrr
         gph300r1covar(:,:,RR)=gph300r1covar(:,:,RR)/gph300var(:,:)
         gph500r1covar(:,:,RR)=gph500r1covar(:,:,RR)/gph500var(:,:)
         gph700r1covar(:,:,RR)=gph700r1covar(:,:,RR)/gph700var(:,:)
         gph850r1covar(:,:,RR)=gph850r1covar(:,:,RR)/gph850var(:,:)
         gph925r1covar(:,:,RR)=gph925r1covar(:,:,RR)/gph925var(:,:)

         t300r1covar(:,:,RR)=t300r1covar(:,:,RR)/t300var(:,:)
         t500r1covar(:,:,RR)=t500r1covar(:,:,RR)/t500var(:,:)
         t700r1covar(:,:,RR)=t700r1covar(:,:,RR)/t700var(:,:)
         t850r1covar(:,:,RR)=t850r1covar(:,:,RR)/t850var(:,:)
         t925r1covar(:,:,RR)=t925r1covar(:,:,RR)/t925var(:,:)

         u300r1covar(:,:,RR)=u300r1covar(:,:,RR)/u300var(:,:)
         u500r1covar(:,:,RR)=u500r1covar(:,:,RR)/u500var(:,:)
         u700r1covar(:,:,RR)=u700r1covar(:,:,RR)/u700var(:,:)
         u850r1covar(:,:,RR)=u850r1covar(:,:,RR)/u850var(:,:)
         u925r1covar(:,:,RR)=u925r1covar(:,:,RR)/u925var(:,:)

         v300r1covar(:,:,RR)=v300r1covar(:,:,RR)/v300var(:,:)
         v500r1covar(:,:,RR)=v500r1covar(:,:,RR)/v500var(:,:)
         v700r1covar(:,:,RR)=v700r1covar(:,:,RR)/v700var(:,:)
         v850r1covar(:,:,RR)=v850r1covar(:,:,RR)/v850var(:,:)
         v925r1covar(:,:,RR)=v925r1covar(:,:,RR)/v925var(:,:)

         q850r1covar(:,:,RR)=q850r1covar(:,:,RR)/q850var(:,:)

         slpr1covar(:,:,RR)=slpr1covar(:,:,RR)/slpvar(:,:)
         q2r1covar(:,:,RR)=q2r1covar(:,:,RR)/q2var(:,:)
         t2r1covar(:,:,RR)=t2r1covar(:,:,RR)/t2var(:,:)
         td2r1covar(:,:,RR)=td2r1covar(:,:,RR)/td2var(:,:)
         u10r1covar(:,:,RR)=u10r1covar(:,:,RR)/u10var(:,:)
         v10r1covar(:,:,RR)=v10r1covar(:,:,RR)/v10var(:,:)

      enddo


C$$$ Fix Units

      q850r1covar(:,:,:)=q850r1covar(:,:,:)/1000.0
      q2r1covar(:,:,:)=q2r1covar(:,:,:)/1000.0
      slpr1covar(:,:,RR)=slpr1covar(:,:,RR)*100.0

      slpmean=slpmean/100.0
      t2mean=t2mean-273.15
      td2mean=td2mean-273.15
      t300mean=t300mean-273.15
      t500mean=t500mean-273.15
      t700mean=t700mean-273.15
      t850mean=t850mean-273.15
      t925mean=t925mean-273.15
      q2mean=q2mean*1000.0
      q850mean=q850mean*1000.0

C$$$ Zero out underground values

      do i=1,mix-1
         do j=1,mjx-1
            if (uggph700(i,j) .GT. 0) then
               gph700r1covar(i,j,:)=9.0E+9
               t700r1covar(i,j,:)=9.0E+9
               u700r1covar(i,j,:)=9.0E+9
               v700r1covar(i,j,:)=9.0E+9
               gph700mean(i,j)=9.0E+9
               t700mean(i,j)=9.0E+9
               gph700targ(i,j,:)=9.0E+9
               t700targ(i,j,:)=9.0E+9
               u700targ(i,j,:)=9.0E+9
               v700targ(i,j,:)=9.0E+9
            endif
            if (uggph850(i,j) .GT. 0) then
               gph850r1covar(i,j,:)=9.0E+9
               t850r1covar(i,j,:)=9.0E+9
               u850r1covar(i,j,:)=9.0E+9
               v850r1covar(i,j,:)=9.0E+9
               q850r1covar(i,j,:)=9.0E+9
               gph850mean(i,j)=9.0E+9
               t850mean(i,j)=9.0E+9
               q850mean(i,j)=9.0E+9
               gph850targ(i,j,:)=9.0E+9
               t850targ(i,j,:)=9.0E+9
               u850targ(i,j,:)=9.0E+9
               v850targ(i,j,:)=9.0E+9
               q850targ(i,j,:)=9.0E+9
            endif
            if (uggph925(i,j) .GT. 0) then
               gph925r1covar(i,j,:)=9.0E+9
               t925r1covar(i,j,:)=9.0E+9
               u925r1covar(i,j,:)=9.0E+9
               v925r1covar(i,j,:)=9.0E+9
               gph925mean(i,j)=9.0E+9
               t925mean(i,j)=9.0E+9
               gph925targ(i,j,:)=9.0E+9
               t925targ(i,j,:)=9.0E+9
               u925targ(i,j,:)=9.0E+9
               v925targ(i,j,:)=9.0E+9
            endif
         enddo
      enddo

C$$$ Assemble single response into compressed data vars
      print*, rrr
      smat(:,:,1)=gph300r1covar(:,:,rrr)
      smat(:,:,2)=gph500r1covar(:,:,rrr)
      smat(:,:,3)=gph700r1covar(:,:,rrr)
      smat(:,:,4)=gph850r1covar(:,:,rrr)
      smat(:,:,5)=t300r1covar(:,:,rrr)
      smat(:,:,6)=t500r1covar(:,:,rrr)
      smat(:,:,7)=t700r1covar(:,:,rrr)
      smat(:,:,8)=t850r1covar(:,:,rrr)
      smat(:,:,9)=t925r1covar(:,:,rrr)
      smat(:,:,10)=u300r1covar(:,:,rrr)
      smat(:,:,11)=u500r1covar(:,:,rrr)
      smat(:,:,12)=u700r1covar(:,:,rrr)
      smat(:,:,13)=u850r1covar(:,:,rrr)
      smat(:,:,14)=u925r1covar(:,:,rrr)
      smat(:,:,15)=v300r1covar(:,:,rrr)
      smat(:,:,16)=v500r1covar(:,:,rrr)
      smat(:,:,17)=v700r1covar(:,:,rrr)
      smat(:,:,18)=v850r1covar(:,:,rrr)
      smat(:,:,19)=v925r1covar(:,:,rrr)
      smat(:,:,20)=q850r1covar(:,:,rrr)
      smat(:,:,21)=slpr1covar(:,:,rrr)
      smat(:,:,22)=t2r1covar(:,:,rrr)
      smat(:,:,23)=q2r1covar(:,:,rrr)
      smat(:,:,24)=td2r1covar(:,:,rrr)
      smat(:,:,25)=u10r1covar(:,:,rrr)
      smat(:,:,26)=v10r1covar(:,:,rrr)

      mmat(:,:,1)=gph300mean(:,:)
      mmat(:,:,2)=gph500mean(:,:)
      mmat(:,:,3)=gph700mean(:,:)
      mmat(:,:,4)=gph850mean(:,:)
      mmat(:,:,5)=t300mean(:,:)
      mmat(:,:,6)=t500mean(:,:)
      mmat(:,:,7)=t700mean(:,:)
      mmat(:,:,8)=t850mean(:,:)
      mmat(:,:,9)=t925mean(:,:)
      mmat(:,:,10)=u300mean(:,:)
      mmat(:,:,11)=u500mean(:,:)
      mmat(:,:,12)=u700mean(:,:)
      mmat(:,:,13)=u850mean(:,:)
      mmat(:,:,14)=u925mean(:,:)
      mmat(:,:,15)=v300mean(:,:)
      mmat(:,:,16)=v500mean(:,:)
      mmat(:,:,17)=v700mean(:,:)
      mmat(:,:,18)=v850mean(:,:)
      mmat(:,:,19)=v925mean(:,:)
      mmat(:,:,20)=q850mean(:,:)
      mmat(:,:,21)=slpmean(:,:)
      mmat(:,:,22)=t2mean(:,:)
      mmat(:,:,23)=q2mean(:,:)
      mmat(:,:,24)=td2mean(:,:)
      mmat(:,:,25)=u10mean(:,:)
      mmat(:,:,26)=v10mean(:,:)

C$$$ Section here for writing esens, targ, means

      call open_file(outfile1,permissw,iunit1)
      call write_variable3d(iunit1,'P_HYD',mix-1,
     &     mjx-1,mkx-1,1,smat)
      call write_variable3d(iunit1,'QICE',mix-1,
     &     mjx-1,mkx-1,1,mmat)
      call close_file(iunit1)

c$$ If specified, save member responses to netcdf file

      if (scatter == .true.) then

C$$$ Save off response vectors for each member
        call open_file(filer, nf_write, fidr)

        call write_variable1d(fidr,'DBZ_AVG',ensnum,1,dbzavgvec)
        call write_variable1d(fidr,'DBZ_MAX',ensnum,1,dbzmaxvec)
        call write_variable1d(fidr,'UH_AVG',ensnum,1,uhavgvec)
        call write_variable1d(fidr,'UH_MAX',ensnum,1,uhmaxvec)
        call write_variable1d(fidr,'PCP',ensnum,1,pcpvec)
        call write_variable1d(fidr,'WSPD_AVG',ensnum,1,windavgvec)
        call write_variable1d(fidr,'UH_COV',ensnum,1,uhcovvec)
        call write_variable1d(fidr,'DBZ_COV',ensnum,1,dbzcovvec)

        call close_file(fidr)

      endif


      end
