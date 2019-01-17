      program meancalc

c$$$ to compile, use Makefile (need netcdf.mod, wrf_tools.mod in dir, compile
c$$$ module_netcdf.f and module_wrf_tools.f if the .mod files don't exist)
c$$$ Also, remove meancalc before compiling with Makefile

c$$$ This program uses the modules to read data from WRF files,
c$$$ calculates means, and writes them back to a 
c$$$ new WRF file for plotting with RIP and for use with other
c$$$ programs.  Need to copy a single ens member into SENSmean.out
c$$$ before running this program.

      USE NETCDF
      USE WRF_TOOLS

      integer mix,mjx,mkx,permissr,permissw,cc
      integer iunitwritemean,rcode
      real, pointer::qvapor(:,:,:)
      real, pointer::dbz(:,:,:)
      real, pointer::totpres(:,:,:)
      real, pointer::t(:,:,:)
      real, pointer::ph(:,:,:)
      real, pointer::phb(:,:,:)
      real, pointer::theta(:,:,:)
      real, pointer::temp(:,:,:)
      real, pointer::rhodry(:,:,:)
      real, pointer::uu(:,:,:)
      real, pointer::vv(:,:,:)
      real, pointer::ww(:,:,:)
      real, pointer::gpttot(:,:,:)
      real, pointer::p(:,:,:)
      real, pointer::phyd(:,:,:)
      real, pointer::pb(:,:,:)
      real, pointer::qrain(:,:,:)
      real, pointer::qcloud(:,:,:)
      real, pointer::qgraup(:,:,:)
      real, pointer::qice(:,:,:)
      real, pointer::qsnow(:,:,:)
      real, pointer::surfptotal(:,:)
      real, pointer::mu(:,:)
      real, pointer::mub(:,:)
      real, pointer::sst(:,:)
      real, pointer::rainnc(:,:)
      real, pointer::snownc(:,:)
      real, pointer::hailnc(:,:)
      real, pointer::graupelnc(:,:)
      real, pointer::rainnum(:,:,:)
      real, pointer::icenum(:,:,:)
      real, pointer::q2(:,:)
      real, pointer::t2(:,:)
      real, pointer::th2(:,:)
      real, pointer::psfc(:,:)
      real, pointer::u10(:,:)
      real, pointer::v10(:,:)
      real, pointer::rdnw(:)
      real, pointer::rdn(:)
      real, pointer::dnw(:)
      real, pointer::dn(:)
      real, pointer::znu(:)
      real, pointer::znw(:)

      real, pointer::qvapornext(:,:,:)
      real, pointer::dbznext(:,:,:)
      real, pointer::totpresnext(:,:,:)
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

      integer, pointer::times(:)
      integer, pointer::levels(:)
      character*67, pointer::infile(:)
      character*20 outfilemean,infileinfo
      real basethet,levee,rdgas,rvgas,preffer
      real kap,prefcont,kapdiv
      integer ensnum,numtimes,iunitread
      integer numpoints,i,j,k,m,ensnummax
      integer timeinit,iunitinfo,ecnt
      character*2 time

      print*, 'Variables declared'

c$$$ Initialize constants

      ecnt=42

      print*,ecnt

      outfilemean='SENSmean.out'
      infileinfo='mem1/SENS1_0.out'
      ensnummax=42
      ensnum=ecnt
      iunitread=10
      iunitinfo=22
      iunitwritemean=11
      basethet=300.
      rdgas=287.
      rvgas=461.6
      preffer=100000.
      kap=2./7.
      prefcont=preffer**(-(kap))
      kapdiv=1./(1.-kap)
      permissw=3
      permissr=2

      call open_file(infileinfo,permissr,iunitinfo)
c$$      numtimes = get_dimlen(iunitinfo,'Time')
c$$   Pull 14 time steps such that sens time can be 1-14
      numtimes=14
      rcode = nf_get_att_int(iunitinfo, nf_global,
     &                    'WEST-EAST_GRID_DIMENSION', mix)
      rcode = nf_get_att_int(iunitinfo, nf_global,
     &                    'SOUTH-NORTH_GRID_DIMENSION', mjx)
      rcode = nf_get_att_int(iunitinfo, nf_global,
     &                    'BOTTOM-TOP_GRID_DIMENSION', mkx)
      allocate(znw(mkx))
      allocate(znu(mkx-1))
c$$      call get_variable1d(iunitinfo,'ZNW',mkx,
c$$     &     1,znw)
c$$      call get_variable1d(iunitinfo,'ZNU',mkx-1,
c$$     &     1,znu)
c$$      call close_file(iunitinfo)

      print*, "Constants initialized"

c$$$ Allocate variable arrays

c$$$ Beginning of loop through each time
      do cc=1,numtimes

      if (cc .EQ. 1) then
      allocate(infile(ensnummax))
      allocate(uu(mix,mjx-1,mkx-1))
      allocate(vv(mix-1,mjx,mkx-1))
      allocate(ww(mix-1,mjx-1,mkx-1))
      allocate(ph(mix-1,mjx-1,mkx))
      allocate(phb(mix-1,mjx-1,mkx))
      allocate(t(mix-1,mjx-1,mkx-1))
      allocate(mu(mix-1,mjx-1))
      allocate(mub(mix-1,mjx-1))
      allocate(p(mix-1,mjx-1,mkx-1))
      allocate(dbz(mix-1,mjx-1,mkx-1))
      allocate(dbznext(mix-1,mjx-1,mkx-1))
      allocate(phyd(mix-1,mjx-1,mkx-1))
      allocate(pb(mix-1,mjx-1,mkx-1))
      allocate(rdnw(mkx-1))
      allocate(rdn(mkx-1))
      allocate(dnw(mkx-1))
      allocate(dn(mkx-1))
      allocate(q2(mix-1,mjx-1))
      allocate(t2(mix-1,mjx-1))
      allocate(th2(mix-1,mjx-1))
      allocate(psfc(mix-1,mjx-1))
      allocate(u10(mix-1,mjx-1))
      allocate(v10(mix-1,mjx-1))
      allocate(qvapor(mix-1,mjx-1,mkx-1))
      allocate(qrain(mix-1,mjx-1,mkx-1))      
      allocate(qcloud(mix-1,mjx-1,mkx-1))
      allocate(qice(mix-1,mjx-1,mkx-1))
      allocate(qsnow(mix-1,mjx-1,mkx-1))      
      allocate(qgraup(mix-1,mjx-1,mkx-1))
      allocate(sst(mix-1,mjx-1))
      allocate(rainnc(mix-1,mjx-1))
      allocate(snownc(mix-1,mjx-1))
      allocate(hailnc(mix-1,mjx-1))
      allocate(graupelnc(mix-1,mjx-1))
      allocate(rainnum(mix-1,mjx-1,mkx-1))
      allocate(icenum(mix-1,mjx-1,mkx-1))
      allocate(totpres(mix-1,mjx-1,mkx-1))
      allocate(theta(mix-1,mjx-1,mkx-1))
      allocate(temp(mix-1,mjx-1,mkx-1))
      allocate(gpttot(mix-1,mjx-1,mkx))
      allocate(rhodry(mix-1,mjx-1,mkx-1))
      allocate(surfptotal(mix-1,mjx-1))
      allocate(uunext(mix,mjx-1,mkx-1))
      allocate(vvnext(mix-1,mjx,mkx-1))
      allocate(wwnext(mix-1,mjx-1,mkx-1))
      allocate(phnext(mix-1,mjx-1,mkx))
      allocate(phbnext(mix-1,mjx-1,mkx))
      allocate(tnext(mix-1,mjx-1,mkx-1))
      allocate(munext(mix-1,mjx-1))
      allocate(mubnext(mix-1,mjx-1))
      allocate(pnext(mix-1,mjx-1,mkx-1))
      allocate(phydnext(mix-1,mjx-1,mkx-1))
      allocate(pbnext(mix-1,mjx-1,mkx-1))
      allocate(rdnwnext(mkx-1))
      allocate(rdnnext(mkx-1))
      allocate(dnwnext(mkx-1))
      allocate(dnnext(mkx-1))
      allocate(q2next(mix-1,mjx-1))
      allocate(t2next(mix-1,mjx-1))
      allocate(th2next(mix-1,mjx-1))
      allocate(psfcnext(mix-1,mjx-1))
      allocate(u10next(mix-1,mjx-1))
      allocate(v10next(mix-1,mjx-1))
      allocate(qvapornext(mix-1,mjx-1,mkx-1))
      allocate(qrainnext(mix-1,mjx-1,mkx-1))      
      allocate(qcloudnext(mix-1,mjx-1,mkx-1))
      allocate(qicenext(mix-1,mjx-1,mkx-1))
      allocate(qsnownext(mix-1,mjx-1,mkx-1))      
      allocate(qgraupnext(mix-1,mjx-1,mkx-1))
      allocate(sstnext(mix-1,mjx-1))
      allocate(rainncnext(mix-1,mjx-1))
      allocate(snowncnext(mix-1,mjx-1))
      allocate(hailncnext(mix-1,mjx-1))
      allocate(graupelncnext(mix-1,mjx-1))
      allocate(rainnumnext(mix-1,mjx-1,mkx-1))
      allocate(icenumnext(mix-1,mjx-1,mkx-1))
      allocate(totpresnext(mix-1,mjx-1,mkx-1))
      allocate(thetanext(mix-1,mjx-1,mkx-1))
      allocate(tempnext(mix-1,mjx-1,mkx-1))
      allocate(gpttotnext(mix-1,mjx-1,mkx))
      allocate(rhodrynext(mix-1,mjx-1,mkx-1))
      allocate(surfptotalnext(mix-1,mjx-1))
      endif

      print*, 'End memory allocations'

      write(time, '(i2)') cc
      print*, cc
      infile = (/ 'mem1/SENS1_' // trim(adjustl(time)) // '.out',
     &     'mem2/SENS2_' // trim(adjustl(time)) // '.out',
     &     'mem3/SENS3_' // trim(adjustl(time)) // '.out',   
     &     'mem4/SENS4_' // trim(adjustl(time)) //'.out',
     &     'mem5/SENS5_' // trim(adjustl(time)) //'.out',
     &     'mem6/SENS6_' // trim(adjustl(time)) //'.out',
     &     'mem7/SENS7_' // trim(adjustl(time)) //'.out',
     &     'mem8/SENS8_' // trim(adjustl(time)) //'.out',
     &     'mem9/SENS9_' // trim(adjustl(time)) //'.out',
     &     'mem10/SENS10_' // trim(adjustl(time)) // '.out',
     &     'mem11/SENS11_' // trim(adjustl(time)) // '.out',
     &     'mem12/SENS12_' // trim(adjustl(time)) // '.out',
     &     'mem13/SENS13_' // trim(adjustl(time)) // '.out',
     &     'mem14/SENS14_' // trim(adjustl(time)) // '.out',
     &     'mem15/SENS15_'// trim(adjustl(time)) //'.out',
     &     'mem16/SENS16_'// trim(adjustl(time)) //'.out',
     &     'mem17/SENS17_'// trim(adjustl(time)) //'.out',
     &     'mem18/SENS18_'// trim(adjustl(time)) //'.out',
     &     'mem19/SENS19_'// trim(adjustl(time)) //'.out',
     &     'mem20/SENS20_'// trim(adjustl(time)) //'.out',
     &     'mem21/SENS21_'// trim(adjustl(time)) //'.out',
     &     'mem22/SENS22_'// trim(adjustl(time)) //'.out',
     &     'mem23/SENS23_'// trim(adjustl(time)) //'.out',
     &     'mem24/SENS24_'// trim(adjustl(time)) //'.out',
     &     'mem25/SENS25_'// trim(adjustl(time)) //'.out',
     &     'mem26/SENS26_'// trim(adjustl(time)) //'.out',
     &     'mem27/SENS27_'// trim(adjustl(time)) //'.out',
     &     'mem28/SENS28_'// trim(adjustl(time)) //'.out',
     &     'mem29/SENS29_'// trim(adjustl(time)) //'.out',
     &     'mem30/SENS30_'// trim(adjustl(time)) //'.out',
     &     'mem31/SENS31_'// trim(adjustl(time)) //'.out',
     &     'mem32/SENS32_'// trim(adjustl(time)) //'.out',
     &     'mem33/SENS33_'// trim(adjustl(time)) //'.out',
     &     'mem34/SENS34_'// trim(adjustl(time)) //'.out',
     &     'mem35/SENS35_'// trim(adjustl(time)) //'.out',
     &     'mem36/SENS36_'// trim(adjustl(time)) //'.out',
     &     'mem37/SENS37_'// trim(adjustl(time)) //'.out',
     &     'mem38/SENS38_'// trim(adjustl(time)) //'.out',
     &     'mem39/SENS39_'// trim(adjustl(time)) //'.out',
     &     'mem40/SENS40_'// trim(adjustl(time)) //'.out',
     &     'mem41/SENS41_'// trim(adjustl(time)) //'.out',
     &     'mem42/SENS42_'// trim(adjustl(time)) //'.out' /)

      print*, "Memory allocated and arrays declared"


c$$$ Read variables here

      print*, "Begin reading variables"

      call open_file(infile(1),permissr,iunitread)
      call get_variable3d(iunitread,'U',mix,
     &     mjx-1,mkx-1,1,uu)
      call get_variable3d(iunitread,'V',mix-1,
     &     mjx,mkx-1,1,vv)
      call get_variable3d(iunitread,'W',mix-1,
     &     mjx-1,mkx-1,1,ww)
      call get_variable3d(iunitread,'PH',mix-1,
     &     mjx-1,mkx,1,ph)
c$$      call get_variable3d(iunitread,'PHB',mix-1,
c$$     &     mjx-1,mkx,0,phb)
      call get_variable3d(iunitread,'T',mix-1,
     &     mjx-1,mkx-1,1,t)
      call get_variable3d(iunitread,'REFL_10CM',mix-1,
     &     mjx-1,mkx-1,1,dbz)
      call get_variable2d(iunitread,'MU',mix-1,
     &     mjx-1,1,mu)
c$$      call get_variable2d(iunitread,'MUB',mix-1,
c$$     &     mjx-1,0,mub)
      call get_variable3d(iunitread,'P',mix-1,
     &     mjx-1,mkx-1,1,p)
c$$      call get_variable3d(iunitread,'P_HYD',mix-1,
c$$     &     mjx-1,mkx-1,0,phyd)
c$$      call get_variable3d(iunitread,'PB',mix-1,
c$$     &     mjx-1,mkx-1,0,pb)
      call get_variable2d(iunitread,'Q2',mix-1,
     &     mjx-1,1,q2)
      call get_variable2d(iunitread,'T2',mix-1,
     &     mjx-1,1,t2)
      call get_variable2d(iunitread,'TH2',mix-1,
     &     mjx-1,1,th2)
      call get_variable2d(iunitread,'PSFC',mix-1,
     &     mjx-1,1,psfc)
      call get_variable2d(iunitread,'U10',mix-1,
     &     mjx-1,1,u10)
      call get_variable2d(iunitread,'V10',mix-1,
     &     mjx-1,1,v10)
      call get_variable3d(iunitread,'QVAPOR',mix-1,
     &     mjx-1,mkx-1,1,qvapor)
      call get_variable3d(iunitread,'QRAIN',mix-1,
     &     mjx-1,mkx-1,1,qrain)
      call get_variable3d(iunitread,'QCLOUD',mix-1,
     &     mjx-1,mkx-1,1,qcloud)
      call get_variable3d(iunitread,'QGRAUP',mix-1,
     &     mjx-1,mkx-1,1,qgraup)
      call get_variable3d(iunitread,'QSNOW',mix-1,
     &     mjx-1,mkx-1,1,qsnow)
c$$      call get_variable3d(iunitread,'QICE',mix-1,
c$$     &     mjx-1,mkx-1,1,qice)
c$$      call get_variale2d(iunitread,'SST',mix-1,
c$$     &     mjx-1,0,sst)
      call get_variable2d(iunitread,'RAINNC',mix-1,
     &     mjx-1,1,rainnc)
c$$      call get_variable2d(iunitread,'SNOWNC',mix-1,
c$$     &     mjx-1,0,snownc)
c$$      call get_variable2d(iunitread,'HAILNC',mix-1,
c$$     &     mjx-1,0,hailnc)
c$$      call get_variable2d(iunitread,'GRAUPELNC',mix-1,
c$$     &     mjx-1,0,graupelnc)
c$$      call get_variable3d(iunitread,'QNRAIN',mix-1,
c$$     &     mjx-1,mkx-1,0,rainnum)
c$$      call get_variable3d(iunitread,'QNICE',mix-1,
c$$     &     mjx-1,mkx-1,0,icenum)
      call close_file(iunitread)

      print*,'got all vars for member'
      print*,1

c$$$ Read in all other ens members, summing data for mean
      do m=2,ensnum
      
      print*,infile(m)

      call open_file(infile(m),permissr,iunitread)
      call get_variable3d(iunitread,'U',mix,
     &     mjx-1,mkx-1,1,uunext)
      call get_variable3d(iunitread,'V',mix-1,
     &     mjx,mkx-1,1,vvnext)
      call get_variable3d(iunitread,'W',mix-1,
     &     mjx-1,mkx-1,1,wwnext)
      call get_variable3d(iunitread,'PH',mix-1,
     &     mjx-1,mkx,1,phnext)
c$$      call get_variable3d(iunitread,'PHB',mix-1,
c$$     &     mjx-1,mkx,0,phbnext)
      call get_variable3d(iunitread,'T',mix-1,
     &     mjx-1,mkx-1,1,tnext)
      call get_variable3d(iunitread,'REFL_10CM',mix-1,
     &     mjx-1,mkx-1,1,dbznext)
      call get_variable2d(iunitread,'MU',mix-1,
     &     mjx-1,1,munext)
c$$      call get_variable2d(iunitread,'MUB',mix-1,
c$$     &     mjx-1,0,mubnext)
      call get_variable3d(iunitread,'P',mix-1,
     &     mjx-1,mkx-1,1,pnext)
c$$      call get_variable3d(iunitread,'P_HYD',mix-1,
c$$     &     mjx-1,mkx-1,0,phydnext)
c$$      call get_variable3d(iunitread,'PB',mix-1,
c$$     &     mjx-1,mkx-1,0,pbnext)
c$$      call get_variable1d(iunitread,'RDNW',mkx-1,
c$$     &     1,rdnwnext)
c$$      call get_variable1d(iunitread,'RDN',mkx-1,
c$$     &     1,rdnnext)
      call get_variable2d(iunitread,'Q2',mix-1,
     &     mjx-1,1,q2next)
      call get_variable2d(iunitread,'T2',mix-1,
     &     mjx-1,1,t2next)
      call get_variable2d(iunitread,'TH2',mix-1,
     &     mjx-1,1,th2next)
      call get_variable2d(iunitread,'PSFC',mix-1,
     &     mjx-1,1,psfcnext)
      call get_variable2d(iunitread,'U10',mix-1,
     &     mjx-1,1,u10next)
      call get_variable2d(iunitread,'V10',mix-1,
     &     mjx-1,1,v10next)
      call get_variable3d(iunitread,'QVAPOR',mix-1,
     &     mjx-1,mkx-1,1,qvapornext)
      call get_variable3d(iunitread,'QRAIN',mix-1,
     &     mjx-1,mkx-1,1,qrainnext)
      call get_variable3d(iunitread,'QCLOUD',mix-1,
     &     mjx-1,mkx-1,1,qcloudnext)
      call get_variable3d(iunitread,'QGRAUP',mix-1,
     &     mjx-1,mkx-1,1,qgraupnext)
      call get_variable3d(iunitread,'QSNOW',mix-1,
     &     mjx-1,mkx-1,1,qsnownext)
c$$      call get_variable3d(iunitread,'QICE',mix-1,
c$$     &     mjx-1,mkx-1,1,qicenext)
c$$      call get_variable2d(iunitread,'SST',mix-1,
c$$     &     mjx-1,0,sstnext)
      call get_variable2d(iunitread,'RAINNC',mix-1,
     &     mjx-1,1,rainncnext)
c$$      call get_variable2d(iunitread,'SNOWNC',mix-1,
c$$     &     mjx-1,0,snowncnext)
c$$      call get_variable2d(iunitread,'HAILNC',mix-1,
c$$     &     mjx-1,0,hailncnext)
c$$      call get_variable2d(iunitread,'GRAUPELNC',mix-1,
c$$     &     mjx-1,0,graupelncnext)
c$$      call get_variable3d(iunitread,'QNRAIN',mix-1,
c$$     &     mjx-1,mkx-1,0,rainnumnext)
c$$      call get_variable3d(iunitread,'QNICE',mix-1,
c$$     &     mjx-1,mkx-1,0,icenumnext)

      call close_file(iunitread)      
      print*,'got all vars for member'
      print*,m

      uu=uu+uunext
      vv=vv+vvnext
      ww=ww+wwnext
      ph=ph+phnext
c$$      phb=phb+phbnext
      t=t+tnext
      dbz=dbz+dbznext
      mu=mu+munext
c$$      mub=mub+mubnext
      p=p+pnext
c$$      phyd=phyd+phydnext
c$$      pb=pb+pbnext
      q2=q2+q2next
      t2=t2+t2next
      th2=th2+th2next
      psfc=psfc+psfcnext
      u10=u10+u10next
      v10=v10+v10next
      qvapor=qvapor+qvapornext
      qrain=qrain+qrainnext
      qcloud=qcloud+qcloudnext
      qgraup=qgraup+qgraupnext
      qsnow=qsnow+qsnownext
c$$      qice=qice+qicenext
c$$      sst=sst+sstnext
      rainnc=rainnc+rainncnext
c$$      snownc=snownc+snowncnext
c$$      rainnum=rainnum+rainnumnext
c$$      icenum=icenum+icenumnext
c$$      hailnc=hailnc+hailncnext
c$$      graupelnc=graupelnc+graupelncnext

c$$$ End of loop through ens members
      enddo

      uu=uu/ensnum
      vv=vv/ensnum
      ww=ww/ensnum
      ph=ph/ensnum
c$$      phb=phb/ensnum
      t=t/ensnum
      mu=mu/ensnum
c$$      mub=mub/ensnum
      p=p/ensnum
      dbz=dbz/ensnum
c$$      phyd=phyd/ensnum
c$$      pb=pb/ensnum
      q2=q2/ensnum
      t2=t2/ensnum
      th2=th2/ensnum
      psfc=psfc/ensnum
      u10=u10/ensnum
      v10=v10/ensnum
      qvapor=qvapor/ensnum
      qrain=qrain/ensnum
      qcloud=qcloud/ensnum
      qsnow=qsnow/ensnum
      qgraup=qgraup/ensnum
c$$      qice=qice/ensnum
c$$      sst=sst/ensnum
      rainnc=rainnc/ensnum
c$$      hailnc=hailnc/ensnum
c$$      snownc=snownc/ensnum
c$$      graupelnc=graupelnc/ensnum
c$$      rainnum=rainnum/ensnum
c$$      icenum=icenum/ensnum

      call open_file(outfilemean,permissw,iunitwritemean)

      call write_variable3d(iunitwritemean,'U',mix,
     &     mjx-1,mkx-1,cc,uu)
      call write_variable3d(iunitwritemean,'V',mix-1,
     &     mjx,mkx-1,cc,vv)
      call write_variable3d(iunitwritemean,'W',mix-1,
     &     mjx-1,mkx-1,cc,ww)
      call write_variable3d(iunitwritemean,'PH',mix-1,
     &     mjx-1,mkx,cc,ph)
c$$      call write_variable3d(iunitwritemean,'PHB',mix-1,
c$$     &     mjx-1,mkx,cc,phb)
      call write_variable3d(iunitwritemean,'T',mix-1,
     &     mjx-1,mkx-1,cc,t)
      call write_variable3d(iunitwritemean,'REFL_10CM',mix-1,
     &     mjx-1,mkx-1,cc,dbz)
      call write_variable2d(iunitwritemean,'MU',mix-1,
     &     mjx-1,cc,mu)
c$$      call write_variable2d(iunitwritemean,'MUB',mix-1,
c$$     &     mjx-1,cc,mub)
      call write_variable3d(iunitwritemean,'P',mix-1,
     &     mjx-1,mkx-1,cc,p)
c$$      call write_variable3d(iunitwritemean,'P_HYD',mix-1,
c$$     &     mjx-1,mkx-1,cc,phyd)
c$$      call write_variable3d(iunitwritemean,'PB',mix-1,
c$$     &     mjx-1,mkx-1,cc,pb)
      call write_variable2d(iunitwritemean,'Q2',mix-1,
     &     mjx-1,cc,q2)
      call write_variable2d(iunitwritemean,'T2',mix-1,
     &     mjx-1,cc,t2)
      call write_variable2d(iunitwritemean,'TH2',mix-1,
     &     mjx-1,cc,th2)
      call write_variable2d(iunitwritemean,'PSFC',mix-1,
     &     mjx-1,cc,psfc)
      call write_variable2d(iunitwritemean,'U10',mix-1,
     &     mjx-1,cc,u10)
      call write_variable2d(iunitwritemean,'V10',mix-1,
     &     mjx-1,cc,v10)
      call write_variable3d(iunitwritemean,'QVAPOR',mix-1,
     &     mjx-1,mkx-1,cc,qvapor)
      call write_variable3d(iunitwritemean,'QRAIN',mix-1,
     &     mjx-1,mkx-1,cc,qrain)
      call write_variable3d(iunitwritemean,'QCLOUD',mix-1,
     &     mjx-1,mkx-1,cc,qcloud)
      call write_variable3d(iunitwritemean,'QSNOW',mix-1,
     &     mjx-1,mkx-1,cc,qsnow)
      call write_variable3d(iunitwritemean,'QGRAUP',mix-1,
     &     mjx-1,mkx-1,cc,qgraup)
c$$      call write_variable3d(iunitwritemean,'QICE',mix-1,
c$$     &     mjx-1,mkx-1,cc,qice)
c$$      call write_variable2d(iunitwritemean,'SST',mix-1,
c$$     &     mjx-1,cc,sst)
      call write_variable2d(iunitwritemean,'RAINNC',mix-1,
     &     mjx-1,cc,rainnc)
c$$      call write_variable2d(iunitwritemean,'SNOWNC',mix-1,
c$$     &     mjx-1,cc,snownc)
c$$      call write_variable2d(iunitwritemean,'HAILNC',mix-1,
c$$     &     mjx-1,cc,hailnc)
c$$      call write_variable2d(iunitwritemean,'GRAUPELNC',mix-1,
c$$     &     mjx-1,cc,graupelnc)
c$$      call write_variable3d(iunitwritemean,'QNRAIN',mix-1,
c$$     &     mjx-1,mkx-1,cc,rainnum)
c$$      call write_variable3d(iunitwritemean,'QNICE',mix-1,
c$$     &     mjx-1,mkx-1,cc,icenum)

      call close_file(iunitwritemean)

c$$$ End of loop through times
      enddo

      print*, "Done"

      end


