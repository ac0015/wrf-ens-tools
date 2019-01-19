C$$$ This program writes out netcdf files of sensitivity and
C$$$ response variables for use in other Python programs.
C$$$ Input files are an ensemble of WRF files in the format
C$$$ as seen below, which follows the naming conventions 
C$$$ laid out in 'rename_wrfout.sh'. Metadata used for 'projt'
C$$$ comes from reference file 'wrfoutREF'. Must provide
C$$$ input text file containing number of ensemble members
C$$$ and time to pull sensitivity variables.

      program sensvector

C$$$ compile with make maecalc (invokes Makefile)

      use wrf_tools
      use netcdf
      use map_utils

      implicit none
 
      integer m
      integer iunitin
      type(proj_info) :: projt
      character*40 modfile,filec,assimfile,qfile,ppfile
      character*40 mods(7)
      real, pointer :: lats(:)
      real, pointer :: lons(:)
      real, pointer :: hgt(:,:)
      real, pointer :: dat(:,:)
      real, pointer :: uu(:,:,:)
      real, pointer :: vv(:,:,:)
      real, pointer::surfptotal(:,:)
      real, pointer::qvapor(:,:,:)
      real, pointer::totpres(:,:,:)
      real, pointer::t(:,:,:)
      real, pointer::ph(:,:,:)
      real, pointer::phb(:,:,:)
      real, pointer::theta(:,:,:)
      real, pointer::temp(:,:,:)
      real, pointer::rhodry(:,:,:)
      real, pointer::gpttot(:,:,:)
      real, pointer::gpheight(:,:,:)
      real, pointer::gpheighthalf(:,:,:)
      real, pointer::p(:,:,:)
      real, pointer::pb(:,:,:)
      real, pointer::mu(:,:)
      real, pointer::mub(:,:)
      real, pointer::psfc(:,:)
      real, pointer::q2(:,:)
      real, pointer :: gph300(:,:)
      real, pointer :: gph500(:,:)  
      real, pointer :: gph700(:,:)  
      real, pointer :: gph850(:,:)  
      real, pointer :: t300(:,:)  
      real, pointer :: t500(:,:)  
      real, pointer :: t700(:,:) 
      real, pointer :: t850(:,:) 
      real, pointer :: t925(:,:) 
      real, pointer :: u300(:,:)
      real, pointer :: u500(:,:)
      real, pointer :: u700(:,:)
      real, pointer :: u850(:,:) 
      real, pointer :: u925(:,:) 
      real, pointer :: v300(:,:)
      real, pointer :: v500(:,:)  
      real, pointer :: v700(:,:)  
      real, pointer :: v850(:,:) 
      real, pointer :: v925(:,:)
      real, pointer :: q850(:,:)
      real, pointer :: slp(:,:)   
      real, pointer :: t2(:,:)  
      real, pointer :: td2(:,:) 
      real, pointer :: u10(:,:) 
      real, pointer :: v10(:,:)
      real, pointer :: etah(:)
      real, pointer :: etaf(:)
      integer rnumqc
      real fbig,uob,vob,grav
      character*5, pointer :: flgs(:,:)
      character*8 goo
      character*2 goo1
      real, pointer :: sensvec(:,:,:,:)
      real wrmsall,trmsall,armsall,wrmsclose,trmsclose
      real wmaeall,tmaeall,wmeall,tmeall
      real armsclose,xd,yd,xdm,ydm
      integer flag,bhflag,lhflag,eotflag,dim,ierr,chkclosett
      integer bhi(50,20),latdomi,mixmm5,mjxmm5,chkcloseww
      real bhr(20,20),time,intval1,intval2,utrue1,vtrue1
      real dir1,spd1,dir2,spd2
      integer e(4),fidc,rnum,idgph300,idgph500,idgph700,idgph850
      integer idt300,idt500,idt700,idt850,idt925
      integer idu300,idu500,idu700,idu850,idu925
      integer idv300,idv500,idv700,idv850,idv925
      integer idq850,idslp,idt2,idq2,idtd2,idu10,idv10
      integer altcount,ppcount,ttcount,uucount,vvcount,aacnt
      integer wcntall,tcntall,acntall,wcntclose,tcntclose
      integer acntclose,iunitassim,rnumass,rnumss,ttcnt,wwcnt
      integer permissr,permissw,iunitmod,mix,mjx,mkx,iunitres
      integer ii,jj,rcode,indom,exnum,iunitout,k,i,j,xdim,ydim
      integer iunitoutu,iunitoutv,iunitoutlat,iunitoutlon
      real modelev,ir,jr,tfit,tfcount,ihr,jhr
      real iclose,jclose,missing2,fmiss
      real ufit,ufcount,vfit,vfcount,val

      integer, pointer :: dims(:)
      integer, pointer :: obstatus(:)
      character*2 senstime
      integer timo, ensnum, timesens, mem, iunitref, ensnummax
      real basethet,levee,rdgas,rvgas,preffer
      real kap,prefcont,kapdiv,lapse,cutclose
      character*60, pointer :: infileSENS(:)
      real, pointer :: uucross(:,:,:), vvcross(:,:,:)

      read*, ensnum
      read*, timesens
 
      iunitres  = 9
      bhflag   = 0
      lhflag   = 1
      eotflag  = 2
      missing2=-9999.0
      fmiss=-9999.0
      fbig=1000.0
      timo=1
      ensnummax=42

      basethet=300.
      rdgas=287.
      rvgas=461.6
      preffer=100000.
      kap=2./7.
      prefcont=preffer**(-(kap))
      kapdiv=1./(1.-kap)
      lapse=0.0065
      grav=9.81
      permissr=2
      permissw=3
      fidc=88
      iunitref=22
      iunitmod=33
      iunitout=44
      ppfile='wrfoutREF'
   
      call set_domain_proj(ppfile, projt)
      call open_file(ppfile,permissr,iunitref)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'WEST-EAST_GRID_DIMENSION', mix)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'SOUTH-NORTH_GRID_DIMENSION', mjx)
      rcode = nf_get_att_int(iunitref, nf_global,
     &                    'BOTTOM-TOP_GRID_DIMENSION', mkx)

      allocate(infileSENS(ensnummax))
      allocate(dims(3))
      allocate(sensvec(mix-1,mjx-1,26,ensnum))
      allocate(etah(mkx-1))
      allocate(etaf(mkx-1))
      allocate(gph300(mix-1,mjx-1))
      allocate(gph500(mix-1,mjx-1))  
      allocate(gph700(mix-1,mjx-1)) 
      allocate(gph850(mix-1,mjx-1))  
      allocate(t300(mix-1,mjx-1)) 
      allocate(t500(mix-1,mjx-1))
      allocate(t700(mix-1,mjx-1))
      allocate(t850(mix-1,mjx-1))
      allocate(t925(mix-1,mjx-1))
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
      allocate(q850(mix-1,mjx-1))
      allocate(slp(mix-1,mjx-1))   
      allocate(t2(mix-1,mjx-1))
      allocate(td2(mix-1,mjx-1))
      allocate(u10(mix-1,mjx-1)) 
      allocate(v10(mix-1,mjx-1))     
      allocate(dat(2,2))
      allocate(ph(mix-1,mjx-1,mkx))
      allocate(phb(mix-1,mjx-1,mkx))
      allocate(t(mix-1,mjx-1,mkx-1))
      allocate(mu(mix-1,mjx-1))
      allocate(mub(mix-1,mjx-1))
      allocate(q2(mix-1,mjx-1))
      allocate(psfc(mix-1,mjx-1))
      allocate(p(mix-1,mjx-1,mkx-1))
      allocate(pb(mix-1,mjx-1,mkx-1))
      allocate(qvapor(mix-1,mjx-1,mkx-1))
      allocate(totpres(mix-1,mjx-1,mkx-1))
      allocate(theta(mix-1,mjx-1,mkx-1))
      allocate(temp(mix-1,mjx-1,mkx-1))
      allocate(gpttot(mix-1,mjx-1,mkx))
      allocate(rhodry(mix-1,mjx-1,mkx-1))
      allocate(surfptotal(mix-1,mjx-1))
      allocate(gpheight(mix-1,mjx-1,mkx))
      allocate(gpheighthalf(mix-1,mjx-1,mkx-1))
      allocate(hgt(mix-1,mjx-1))
      allocate(uu(mix-1,mjx-1,mkx-1))
      allocate(vv(mix-1,mjx-1,mkx-1))
      allocate(uucross(mix-1,mjx-1,mkx-1))
      allocate(vvcross(mix-1,mjx-1,mkx-1))

      call get_variable1d(iunitref,'ZNW',mkx,
     &     1,etaf)
      call get_variable1d(iunitref,'ZNU',mkx-1,
     &     1,etah)
      call get_variable2d(iunitref,'HGT',mix-1,
     &     mjx-1,1,hgt)
      call get_variable3d(iunitref,'PHB',mix-1,
     &     mjx-1,mkx,timo,phb)
      call get_variable2d(iunitref,'MUB',mix-1,
     &     mjx-1,timo,mub)
      call close_file(iunitref)

      write(senstime, '(i2)') timesens
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

C$$$ Write netcdf file of sensitivity variables for user-specified time

        filec = 'SENSvals.nc'

        rcode = nf_create(filec, nf_clobber, fidc)

        print*, "Created file ", filec

        rcode = nf_def_dim(fidc, 'xdim', mix-1, xdim)
        rcode = nf_def_dim(fidc, 'ydim', mjx-1, ydim)
        rcode = nf_def_dim(fidc, 'member', ensnum, mem)        

        dims(1) = xdim
        dims(2) = ydim
        dims(3) = mem

        rcode = nf_def_var(fidc,'GPH_300',nf_float,
     &          3,dims,idgph300)
        rcode = nf_put_att_text(fidc, idgph300, 'description', 33,
     &                 '300 hPa geopotential heights at senstime')
        rcode = nf_put_att_text(fidc, idgph300, 'units', 9, 'meters')
        rcode = nf_put_att_real(fidc,idgph300,'_FillValue',nf_float,1,
     &                     fmiss)       

        rcode = nf_def_var(fidc,'GPH_500',nf_float,
     &          3,dims,idgph500)
        rcode = nf_put_att_text(fidc, idgph500, 'description', 33,
     &                 '500 hPa geopotential heights at senstime')
        rcode = nf_put_att_text(fidc, idgph500, 'units', 9, 'meters')
        rcode = nf_put_att_real(fidc,idgph500,'_FillValue',nf_float,1,
     &                     fmiss)  

        rcode = nf_def_var(fidc,'GPH_700',nf_float,
     &          3,dims,idgph700)
        rcode = nf_put_att_text(fidc, idgph700, 'description', 33,
     &                 '700 hPa geopotential heights at senstime')
        rcode = nf_put_att_text(fidc, idgph700, 'units', 9, 'meters')
        rcode = nf_put_att_real(fidc,idgph700,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'GPH_850',nf_float,
     &          3,dims,idgph850)
        rcode = nf_put_att_text(fidc, idgph850, 'description', 33,
     &                 '850 hPa geopotential heights at senstime')
        rcode = nf_put_att_text(fidc, idgph850, 'units', 9, 'meters')
        rcode = nf_put_att_real(fidc,idgph850,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'T_300',nf_float,
     &          3,dims,idt300)
        rcode = nf_put_att_text(fidc, idt300, 'description', 33,
     &                 '300 hPa Temperature')
        rcode = nf_put_att_text(fidc, idt300, 'units', 9, 'Kelvin')
        rcode = nf_put_att_real(fidc,idt300,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'T_500',nf_float,
     &          3,dims,idt500)
        rcode = nf_put_att_text(fidc, idt500, 'description', 33,
     &                 '500 hPa Temperature')
        rcode = nf_put_att_text(fidc, idt500, 'units', 9, 'Kelvin')
        rcode = nf_put_att_real(fidc,idt500,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'T_700',nf_float,
     &          3,dims,idt700)
        rcode = nf_put_att_text(fidc, idt700, 'description', 33,
     &                 '700 hPa Temperature')
        rcode = nf_put_att_text(fidc, idt700, 'units', 9, 'Kelvin')
        rcode = nf_put_att_real(fidc,idt700,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'T_850',nf_float,
     &          3,dims,idt850)
        rcode = nf_put_att_text(fidc, idt850, 'description', 33,
     &                 '850 hPa Temperature')
        rcode = nf_put_att_text(fidc, idt850, 'units', 9, 'Kelvin')
        rcode = nf_put_att_real(fidc,idt850,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'T_925',nf_float,
     &          3,dims,idt925)
        rcode = nf_put_att_text(fidc, idt925, 'description', 33,
     &                 '925 hPa Temperature')
        rcode = nf_put_att_text(fidc, idt925, 'units', 9, 'Kelvin')
        rcode = nf_put_att_real(fidc,idt925,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'U_300',nf_float,
     &          3,dims,idu300)
        rcode = nf_put_att_text(fidc, idu300, 'description', 33,
     &                 '300 hPa U-wind')
        rcode = nf_put_att_text(fidc,idu300,'units', 9, 'meters/second')
        rcode = nf_put_att_real(fidc,idu300,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'U_500',nf_float,
     &          3,dims,idu500)
        rcode = nf_put_att_text(fidc, idu500, 'description', 33,
     &                 '500 hPa U-wind')
        rcode = nf_put_att_text(fidc,idu500, 'units',9, 'meters/second')
        rcode = nf_put_att_real(fidc,idu500,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'U_700',nf_float,
     &          3,dims,idu700)
        rcode = nf_put_att_text(fidc, idu700, 'description', 33,
     &                 '700 hPa U-wind')
        rcode = nf_put_att_text(fidc,idu700,'units', 9, 'meters/second')
        rcode = nf_put_att_real(fidc,idu700,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'U_850',nf_float,
     &          3,dims,idu850)
        rcode = nf_put_att_text(fidc, idu850, 'description', 33,
     &                 '850 hPa U-wind')
        rcode = nf_put_att_text(fidc,idu850, 'units', 9,'meters/second')
        rcode = nf_put_att_real(fidc,13,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'U_925',nf_float,
     &          3,dims,idu925)
        rcode = nf_put_att_text(fidc, idu925, 'description', 33,
     &                 '925 hPa U-wind')
        rcode = nf_put_att_text(fidc,idu925,'units', 9, 'meters/second')
        rcode = nf_put_att_real(fidc,idu925,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'V_300',nf_float,
     &          3,dims,idv300)
        rcode = nf_put_att_text(fidc, idv300, 'description', 33,
     &                 '300 hPa V-wind')
        rcode = nf_put_att_text(fidc,idv300,'units', 9, 'meters/second')
        rcode = nf_put_att_real(fidc,idv300,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'V_500',nf_float,
     &          3,dims,idv500)
        rcode = nf_put_att_text(fidc, idv500, 'description', 33,
     &                 '500 hPa V-wind')
        rcode = nf_put_att_text(fidc,idv500,'units', 9, 'meters/second')
        rcode = nf_put_att_real(fidc,idv500,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'V_700',nf_float,
     &          3,dims,idv700)
        rcode = nf_put_att_text(fidc, idv700, 'description', 33,
     &                 '700 hPa V-wind')
        rcode = nf_put_att_text(fidc,idv700,'units', 9, 'meters/second')
        rcode = nf_put_att_real(fidc,idv700,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'V_850',nf_float,
     &          3,dims,idv850)
        rcode = nf_put_att_text(fidc,idv850, 'description', 33,
     &                 '850 hPa V-wind')
        rcode = nf_put_att_text(fidc,idv850,'units',19,'meters/second')
        rcode = nf_put_att_real(fidc,idv850,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'V_925',nf_float,
     &          3,dims,idv925)
        rcode = nf_put_att_text(fidc, idv925, 'description', 33,
     &                 '925 hPa V-wind')
        rcode = nf_put_att_text(fidc,idv925,'units',19,'meters/second')
        rcode = nf_put_att_real(fidc,idv925,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'Q_850',nf_float,
     &          3,dims,idq850)
        rcode = nf_put_att_text(fidc, idq850, 'description', 33,
     &                 '850 hPa Mixing Ratio')
        rcode = nf_put_att_text(fidc, idq850, 'units', 9, 'kg / kg')
        rcode = nf_put_att_real(fidc,idq850,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'SLP',nf_float,
     &          3,dims,idslp)
        rcode = nf_put_att_text(fidc, idslp, 'description', 33,
     &                 'Sea Level Pressure')
        rcode = nf_put_att_text(fidc, idslp, 'units', 9, 'Pascals')
        rcode = nf_put_att_real(fidc,idslp,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'T2',nf_float,
     &          3,dims,idt2)
        rcode = nf_put_att_text(fidc, idt2, 'description', 33,
     &                 '2 meter Temperature')
        rcode = nf_put_att_text(fidc, idt2, 'units', 15,'Kelvin')
        rcode = nf_put_att_real(fidc,idt2,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'Q2',nf_float,
     &          3,dims,idq2)
        rcode = nf_put_att_text(fidc, idq2, 'description', 33,
     &                 '2 meter Mixing Ratio')
        rcode = nf_put_att_text(fidc, idq2, 'units', 9, 'kg / kg')
        rcode = nf_put_att_real(fidc,idq2,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'TD2',nf_float,
     &          3,dims,idtd2)
        rcode = nf_put_att_text(fidc, idtd2, 'description', 33,
     &                 '2 meter Dewpoint')
        rcode = nf_put_att_text(fidc, idtd2, 'units', 9, 'Fahrenheit')
        rcode = nf_put_att_real(fidc,idtd2,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'U10',nf_float,
     &          3,dims,idu10)
        rcode = nf_put_att_text(fidc, idu10, 'description', 33,
     &                 '10 meter U-wind')
        rcode = nf_put_att_text(fidc, idu10, 'units', 9,'meters/second')
        rcode = nf_put_att_real(fidc,idu10,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'V10',nf_float,
     &          3,dims,idv10)
        rcode = nf_put_att_text(fidc, idv10, 'description', 33,
     &                 '10 meter V-wind')
        rcode = nf_put_att_text(fidc, idv10, 'units', 9,'meters/second')
        rcode = nf_put_att_real(fidc,idv10,'_FillValue',nf_float,1,
     &                     fmiss)

      call close_file(fidc)


      do m=1, ensnum
         print*, "Current ensemble member file path: ", infileSENS(m)
         call open_file(infileSENS(m),permissr,iunitmod)
         call get_variable3d(iunitmod,'PH',mix-1,
     &     mjx-1,mkx,timo,ph)
         call get_variable3d(iunitmod,'T',mix-1,
     &     mjx-1,mkx-1,timo,t)
         call get_variable3d(iunitmod,'U',mix-1,
     &     mjx-1,mkx-1,timo,uu)
         call get_variable3d(iunitmod,'V',mix-1,
     &     mjx-1,mkx-1,timo,vv)
         call get_variable3d(iunitmod,'QVAPOR',mix-1,
     &     mjx-1,mkx-1,timo,qvapor)
         call get_variable2d(iunitmod,'MU',mix-1,
     &     mjx-1,timo,mu)
         call get_variable2d(iunitmod,'U10',mix-1,
     &     mjx-1,timo,u10)
         call get_variable2d(iunitmod,'V10',mix-1,
     &     mjx-1,timo,v10)
         call get_variable2d(iunitmod,'T2',mix-1,
     &     mjx-1,timo,t2)
         call get_variable2d(iunitmod,'Q2',mix-1,
     &     mjx-1,timo,q2)
         call close_file(iunitmod)

      print*, "Begin post-processing"

c$$$ Get pressure, temp, gph on half levels

c$$$  Destagger winds

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

c$$$ Get geopotential height on half levels

         do i=1,mix-1
            do j=1,mjx-1
               call destag_zstag(etah, etaf, mkx-1, 
     &           gpheight(i,j,:),gpheighthalf(i,j,:))
            enddo
         enddo

c$$$  Next calculate dry total surface pressure

         do i=1,mix-1
            do j=1,mjx-1
               surfptotal(i,j)=mu(i,j)+
     &                      mub(i,j)
            enddo
         enddo

      do k=1,mkx-1
         levee=etaf(k+1)-etaf(k)
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
     &            rhodry(i,j,k)*rdgas*(1.+(rvgas/rdgas)*
     &            qvapor(i,j,k)))**kapdiv
               enddo
            enddo
         enddo

c$$$ Now the temperature

         do k=1,mkx-1
            do j=1,mjx-1
               do i=1,mix-1
                  temp(i,j,k)=theta(i,j,k)*
     &            ((preffer/totpres(i,j,k))**-kap)
               enddo
            enddo
         enddo

C$$$ Calculate surface pressure

         do i=1,mix-1
            do j=1,mjx-1
               psfc(i,j)=totpres(i,j,1)*
     &              (2.71828**((gpheighthalf(i,j,1)
     &              -hgt(i,j))/(29.3*temp(i,j,1)
     &              *(1+qvapor(i,j,1)))))
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
            enddo
         enddo

         print*,gph300(30,50)

c$$$ Now calculate temperature from theta,totpres

         do k=1,mkx-1
            do j=1,mjx-1
               do i=1,mix-1
                  temp(i,j,k)=theta(i,j,k)*
     &            ((preffer/totpres(i,j,k))**-kap)
               enddo
            enddo
         enddo

         do i=1,mix-1
            do j=1,mjx-1
               psfc(i,j)=totpres(i,j,1)*
     &              (2.71828**((gpheighthalf(i,j,1)
     &              -hgt(i,j))/(29.3*temp(i,j,1)
     &              *(1+qvapor(i,j,1)))))
            enddo
         enddo

c$$$ Calculate SLP

         do i=1,mix-1
            do j=1,mjx-1
               slp(i,j)=slp_standard_atmos(t2(i,j),
     &              psfc(i,j),q2(i,j),hgt(i,j))
            enddo
         enddo

c$$$ Calculate dew point

         do i=1,mix-1
            do j=1,mjx-1
               td2(i,j)=mixrat_to_tdew(q2(i,j),
     &              psfc(i,j))
            enddo
         enddo

c$$$ Interpolate U, V, T to desired pressure levels

         do i=1,mix-1
            do j=1,mjx-1
               t300(i,j)=interp_pres(temp(i,j,:),
     &              totpres(i,j,:),30000.,mkx-1)
               t500(i,j)=interp_pres(temp(i,j,:),
     &              totpres(i,j,:),50000.,mkx-1)
               t700(i,j)=interp_pres(temp(i,j,:),
     &              totpres(i,j,:),70000.,mkx-1)
               t850(i,j)=interp_pres(temp(i,j,:),
     &              totpres(i,j,:),85000.,mkx-1)
               t925(i,j)=interp_pres(temp(i,j,:),
     &              totpres(i,j,:),92500.,mkx-1)
         enddo
      enddo

         do i=1,mix-1
            do j=1,mjx-1
               u300(i,j)=interp_pres(uucross(i,j,:),
     &              totpres(i,j,:),30000.,mkx-1)
               u500(i,j)=interp_pres(uucross(i,j,:),
     &              totpres(i,j,:),50000.,mkx-1)
               u700(i,j)=interp_pres(uucross(i,j,:),
     &              totpres(i,j,:),70000.,mkx-1)
               u850(i,j)=interp_pres(uucross(i,j,:),
     &              totpres(i,j,:),85000.,mkx-1)
               u925(i,j)=interp_pres(uucross(i,j,:),
     &              totpres(i,j,:),92500.,mkx-1)
            enddo
         enddo

         do i=1,mix-1
            do j=1,mjx-1
               v300(i,j)=interp_pres(vvcross(i,j,:),
     &              totpres(i,j,:),30000.,mkx-1)
               v500(i,j)=interp_pres(vvcross(i,j,:),
     &              totpres(i,j,:),50000.,mkx-1)
               v700(i,j)=interp_pres(vvcross(i,j,:),
     &              totpres(i,j,:),70000.,mkx-1)
               v850(i,j)=interp_pres(vvcross(i,j,:),
     &              totpres(i,j,:),85000.,mkx-1)
               v925(i,j)=interp_pres(vvcross(i,j,:),
     &              totpres(i,j,:),92500.,mkx-1)
            enddo
         enddo

         do i=1,mix-1
            do j=1,mjx-1
               q850(i,j)=interp_pres(qvapor(i,j,:),
     &              totpres(i,j,:),85000.,mkx-1)
            enddo
         enddo


c$$$ Now write out data

         call open_file(filec, nf_write, fidc)

         call write_variable2d(fidc,'GPH_300',mix-1,
     &        mjx-1,m,gph300)
         call write_variable2d(fidc,'GPH_500',mix-1,
     &        mjx-1,m,gph500)
         call write_variable2d(fidc,'GPH_700',mix-1,
     &        mjx-1,m,gph700)
         call write_variable2d(fidc,'GPH_850',mix-1,
     &        mjx-1,m,gph850)
         call write_variable2d(fidc,'T_300',mix-1,
     &        mjx-1,m,t300)
         call write_variable2d(fidc,'T_500',mix-1,
     &        mjx-1,m,t500)
         call write_variable2d(fidc,'T_700',mix-1,
     &        mjx-1,m,t700)
         call write_variable2d(fidc,'T_850',mix-1,
     &        mjx-1,m,t850)
         call write_variable2d(fidc,'T_925',mix-1,
     &        mjx-1,m,t925)
         call write_variable2d(fidc,'U_300',mix-1,
     &        mjx-1,m,u300)
         call write_variable2d(fidc,'U_500',mix-1,
     &        mjx-1,m,u500)
         call write_variable2d(fidc,'U_700',mix-1,
     &        mjx-1,m,u700)
         call write_variable2d(fidc,'U_850',mix-1,
     &        mjx-1,m,u850)
         call write_variable2d(fidc,'U_925',mix-1,
     &        mjx-1,m,u925)
         call write_variable2d(fidc,'V_300',mix-1,
     &        mjx-1,m,v300)
         call write_variable2d(fidc,'V_500',mix-1,
     &        mjx-1,m,v500)
         call write_variable2d(fidc,'V_700',mix-1,
     &        mjx-1,m,v700)
         call write_variable2d(fidc,'V_850',mix-1,
     &        mjx-1,m,v850)
         call write_variable2d(fidc,'V_925',mix-1,
     &        mjx-1,m,v925)
         call write_variable2d(fidc,'Q_850',mix-1,
     &        mjx-1,m,q850)
         call write_variable2d(fidc,'SLP',mix-1,
     &        mjx-1,m,slp)
         call write_variable2d(fidc,'T2',mix-1,
     &        mjx-1,m,t2)
         call write_variable2d(fidc,'Q2',mix-1,
     &        mjx-1,m,q2)
         call write_variable2d(fidc,'TD2',mix-1,
     &        mjx-1,m,td2)
         call write_variable2d(fidc,'U10',mix-1,
     &        mjx-1,m,u10)
         call write_variable2d(fidc,'V10',mix-1,
     &        mjx-1,m,v10)

         call close_file(fidc)
 
c        sensvec(:,:,m,1) = gph300(:,:)
c	 sensvec(:,:,m,2) = gph500(:,:)
c	 sensvec(:,:,m,3) = gph700(:,:)
c	 sensvec(:,:,m,4) = gph850(:,:)
c	 sensvec(:,:,m,5) = t300(:,:)
c	 sensvec(:,:,m,6) = t500(:,:)
c	 sensvec(:,:,m,7) = t700(:,:)
c	 sensvec(:,:,m,8) = t850(:,:)
c	 sensvec(:,:,m,9) = t925(:,:)
c	 sensvec(:,:,m,10) = u300(:,:)
c	 sensvec(:,:,m,11) = u500(:,:)
c	 sensvec(:,:,m,12) = u700(:,:)
c	 sensvec(:,:,m,13) = u850(:,:)
c	 sensvec(:,:,m,14) = u925(:,:)
c	 sensvec(:,:,m,15) = v300(:,:)
c	 sensvec(:,:,m,16) = v500(:,:)
c	 sensvec(:,:,m,17) = v700(:,:)
c	 sensvec(:,:,m,18) = v850(:,:)
c	 sensvec(:,:,m,19) = v925(:,:)
c	 sensvec(:,:,m,20) = q850(:,:)
c	 sensvec(:,:,m,21) = slp(:,:)
c	 sensvec(:,:,m,22) = t2(:,:)
c	 sensvec(:,:,m,23) = q2(:,:)
c	 sensvec(:,:,m,24) = td2(:,:)
c	 sensvec(:,:,m,25) = u10(:,:)
c	 sensvec(:,:,m,26) = v10(:,:)

         print*, "Done with member ", m
        
         enddo
         
       end

