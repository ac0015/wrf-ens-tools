      program probcalc

c$$$ to compile, use Makefile (need netcdf.mod, wrf_tools.mod in dir, compile
c$$$ module_netcdf.f and module_wrf_tools.f if the .mod files don't exist)
c$$$ Also, remove probcalc before compiling with Makefile

c$$$ This program uses the modules to read data from WRF files,
c$$$ calculates probabilities, and writes them back to a 
c$$$ new WRF file for plotting and verification.  
c$$$ Need to copy a single ens member to output file name
c$$$ before running this program. 

c$$$ Takes an input file path from cmd line that contains
c$$$ an integer of the ensemble size, list of members to
c$$$ calculate probabilities for, forecast hour integer
c$$$ for which the probabilities will be valid, a double
c$$$ describing the neighborhood value to use (in km), and
c$$$ finally a string for the outputfile name.

c$$$ Brian Ancell
c$$$ Modifications made by Austin Coleman - 2018
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      USE NETCDF
      USE WRF_TOOLS

      integer mix,mjx,mkx,permissr,permissw,cc
      integer iunitwritemean,rcode,dbcnt
      real, pointer::dbz(:,:,:)
      real, pointer::uphel(:,:)
      real, pointer::wsp(:,:)
      real, pointer::dbzprob(:,:)
      real, pointer::uhprob1(:,:)
      real, pointer::uhprob2(:,:)
      real, pointer::uhprob3(:,:)
      real, pointer::wspprob(:,:)
      real, pointer::probmat(:,:,:)
      character*67, pointer::infile(:,:)
      character*40 outfilemean,infileinfo
      real basethet,levee,rdgas,rvgas,preffer
      real kap,prefcont,kapdiv,distz,rx
      integer ensnum,numtimes,iunitread,uhcnt1,uhcnt2,uhcnt3,wcnt
      integer numpoints,i,j,k,m,t,f,ensnummax,iii,jjj
      integer timeinit,iunitinfo,ecnt,ibeg,iend,jbeg,jend
      integer, pointer::mems(:)
      character*2 fhr,time
      character*2, pointer::sub_mem(:)
      real nbr

c$$$ Initialize constants

      read*,ecnt
      allocate(mems(ecnt))
      read*,mems(:)
      print*,mems
      read*,fhr
c$$$ Read in neighborhood in km
      read*,nbr
c$$$ Read outfilepath name
      read*,outfilemean

      ensnummax=50
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
      cc=1

      call open_file(outfilemean,permissr,iunitinfo)
      rcode = nf_get_att_int(iunitinfo, nf_global,
     &                    'WEST-EAST_GRID_DIMENSION', mix)
      rcode = nf_get_att_int(iunitinfo, nf_global,
     &                    'SOUTH-NORTH_GRID_DIMENSION', mjx)
      rcode = nf_get_att_int(iunitinfo, nf_global,
     &                    'BOTTOM-TOP_GRID_DIMENSION', mkx)
      rcode = nf_get_att_int(iunitinfo, nf_global, 'DX', rx)
      call close_file(iunitinfo)

c$$$ Allocate variable arrays

      allocate(infile(6,ensnum))
      allocate(sub_mem(ensnum))

      allocate(dbz(mix-1,mjx-1,mkx-1))
      allocate(uphel(mix-1,mjx-1))
      allocate(wsp(mix-1,mjx-1))
      allocate(dbzprob(mix-1,mjx-1))
      allocate(uhprob1(mix-1,mjx-1))
      allocate(uhprob2(mix-1,mjx-1))
      allocate(uhprob3(mix-1,mjx-1))
      allocate(wspprob(mix-1,mjx-1))
      allocate(probmat(mix-1,mjx-1,mkx-1))

      print*,'check1'
   
      sub_mem(:) = '0'
      infile(:,:) = ''
      print*, fhr
      read(fhr, '(i2)') f
      print*, f
      do t=f-5,f
      print*, t
      do m=1,ensnum
      write(sub_mem(m), '(i2)') mems(m)
      write(time, '(i2)') t
c$$   get t back into i=1,2,3,4,5,6 etc. form
c$$    by adding six and subtracting final hr
c$$    print*, t-f+6
      infile(t-f+6,m) = ('../mem' // trim(adjustl(sub_mem(m))) //
     & '/R' // trim(adjustl(sub_mem(m))) // '_' //
     & trim(adjustl(time)) // '.out')
      print*, infile(t-f+6, m)
      enddo
      enddo
      print*,'Built sub_mem'  
      print*,infile   

      dbzprob(:,:)=0.0
      uhprob1(:,:)=0.0
      uhprob2(:,:)=0.0
      uhprob3(:,:)=0.0
      wspprob(:,:)=0.0

c$$$ Read in all ens members for probs
      do t=1,6
      do m=1,ensnum
      
      print*, 'Reading member', m
      print*, ' at time', t
      print*, infile(t,m)
      call open_file(infile(t,m),permissr,iunitread)
      call get_variable3d(iunitread,'REFL_10CM',mix-1,
     &     mjx-1,mkx-1,cc,dbz)
      call get_variable2d(iunitread,'UP_HELI_MAX',mix-1,
     &     mjx-1,cc,uphel)
      call get_variable2d(iunitread,'WSPD10MAX',mix-1,
     &     mjx-1,cc,wsp)
      call close_file(iunitread)
      do i=1,mix-1
         do j=1,mjx-1

            if (i .le. 10) then
               ibeg=1
               iend=i+9
            elseif (i .ge. mix-10) then
               ibeg=mix-10
               iend=mix-1
            else
               ibeg=i-9
               iend=i+9
            endif
            if (j .le. 10) then
               jbeg=1
               jend=j+9
            elseif (j .ge. mjx-10) then
               jbeg=mjx-10
               jend=mjx-1
            else
               jbeg=j-9
               jend=j+9
            endif

            dbcnt=0
            uhcnt1=0
            uhcnt2=0
            uhcnt3=0
            wcnt=0
            do iii=ibeg,iend
               do jjj=jbeg,jend
            distz=(((i-iii)*rx)**2 + ((j-jjj)*rx)**2)**0.5

            if (distz.le.nbr) then

            if (dbz(iii,jjj,1) .GE. 40.) then
               dbcnt=dbcnt+1
            endif  

            if (uphel(iii,jjj) .GE. 25.) then
               uhcnt1=uhcnt1+1
            endif      

            if (uphel(iii,jjj) .GE. 40.) then
               uhcnt2=uhcnt2+1
            endif    

            if (uphel(iii,jjj) .GE. 100.) then
                uhcnt3=uhcnt3+1
            endif

            if (wsp(iii,jjj) .GE. 17.8816) then
               wcnt=wcnt+1
            endif    

            endif

               enddo
            enddo

            if (dbcnt.gt.0) then
               dbzprob(i,j)=dbzprob(i,j)+1.0
            endif
            if (uhcnt1.gt.0) then
               uhprob1(i,j)=uhprob1(i,j)+1.0
            endif
            if (uhcnt2.gt.0) then
               uhprob2(i,j)=uhprob2(i,j)+1.0
            endif
            if (uhcnt3.gt.0) then
                uhprob3(i,j)=uhprob3(i,j)+1.0
            endif
            if (wcnt.gt.0) then
               wspprob(i,j)=wspprob(i,j)+1.0
            endif

         enddo
      enddo

      print*,'got all vars for prob for member'
      print*,m

c$$$ End of loop through ens members
      enddo
c$$$ End of loop through times
      enddo
      dbzprob(:,:)=(dbzprob(:,:)/ensnum)*100
      uhprob1(:,:)=(uhprob1(:,:)/ensnum)*100
      uhprob2(:,:)=(uhprob2(:,:)/ensnum)*100
      uhprob3(:,:)=(uhprob3(:,:)/ensnum)*100
      wspprob(:,:)=(wspprob(:,:)/ensnum)*100

      probmat(:,:,:)=0.0
      do i=1,mix-1
         do j=1,mjx-1
            probmat(i,j,1)=dbzprob(i,j)
            probmat(i,j,2)=uhprob1(i,j)
            probmat(i,j,3)=uhprob2(i,j)
            probmat(i,j,4)=uhprob3(i,j)
            probmat(i,j,5)=wspprob(i,j)
         enddo
      enddo

c$$$ If probabilities are 100, set to 99.8 for plotting

      do i=1,mix-1
         do j=1,mjx-1
            do k=1,5
            if (probmat(i,j,k) .GT. 99.9) then
               probmat(i,j,k)=99.8
            endif
            enddo
         enddo
      enddo

      print*, MAXVAL(probmat)

      call open_file(outfilemean,permissw,iunitwritemean)
      call write_variable3d(iunitwritemean,'P_HYD',mix-1,
     &     mjx-1,mkx-1,cc,probmat)
      call close_file(iunitwritemean)


      end


