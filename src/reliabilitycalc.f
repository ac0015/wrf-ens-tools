c$$$ This program calculates probabilistic reliability within a response
c$$$ box (for plotting attribute diagrams) of a forecast given a probability file
c$$$ and a grid of observations as binary hits or misses (0's and 1's).
c$$$ Returns a netCDF file containing the probability bins, forecast
c$$$ frequencies, and observation hit rates.

        program reliabilitycalc

c$$$ Need Makefile to compile

        use wrf_tools
        use netcdf
        use map_utils

        implicit none

        type(proj_info) :: projt
        integer i,j,k,iii,jjj,kkk,mix,mjx,mkx,fhr,fidc
        integer timo,iunitprob,iunitob,ibeg,iend,jbeg,jend
        integer permissr, permissw, nbins, probvarind, bin
        integer rcode, p, idprobbins, idfcstfreq, idobhitrate
        integer bufrsizehint, obid, obcnt, fcsthr
        real thresh, nbrhd, dx, num_nbrhd_pts, fmiss, prob, distz
        character*70 probvar, obvar, sixhr, name
        character*200 probfile, obfile, outfile, variable
        logical sixhour
        real, pointer :: lats(:,:)
        real, pointer :: lons(:,:)
        real, pointer :: fcstprobs(:,:,:)
        real, pointer :: obs(:,:)
        real, pointer :: prob_bins(:)
        real, pointer :: fcst_freq(:)
        real, pointer :: ob_hit_rate(:)
        real, pointer :: rbox_bounds(:)

c$$$ Read input namelist

        allocate(rbox_bounds(4))
        read*, probfile
        print*, "Probfile ", probfile
        read*, obfile
        print*, "Ob file ", obfile
        read*, outfile
        print*, "Outpath ", outfile
        read*, fhr
        print*, "fhr ", fhr
        read*, variable
        print*, "Rel var ", variable
        read*, thresh
        read*, sixhour
        read*, nbrhd
        print*, "Neighborhood (km)", nbrhd
        read*, rbox_bounds(:)

        permissr = 2
        permissw = 3
        iunitprob = 5
        iunitob = 6
        fidc = 88
        timo = 1
        fmiss = 9.0E9
        probvar = "P_HYD"
        obvar = "P_HYD"
        nbins = 10

c$$$ Determine probability index with the chosen prob variable

        if (variable .eq. "updraft_helicity") then
           if (thresh .eq. 25) then
              probvarind = 2
           elseif (thresh .eq. 40) then
              probvarind = 3
           elseif (thresh .eq. 100) then
              probvarind = 4
           endif
        elseif (variable .eq. "reflectivity") then
           if (thresh .eq. 40) then
              probvarind = 1
           elseif (thresh .eq. 50) then
              probvarind = 5
           endif
        endif

        print*, "reliability of", variable, ">", thresh
        print*, probvarind

c$$$ Open probability file for reference information

        call set_domain_proj(probfile, projt)
        call open_file(probfile,permissr,iunitprob)
        rcode = nf_get_att_int(iunitprob, nf_global,
     &                    'WEST-EAST_GRID_DIMENSION', mix)
        rcode = nf_get_att_int(iunitprob, nf_global,
     &                    'SOUTH-NORTH_GRID_DIMENSION', mjx)
        rcode = nf_get_att_int(iunitprob, nf_global,
     &                    'BOTTOM-TOP_GRID_DIMENSION', mkx)
        rcode = nf_get_att_real(iunitprob, nf_global,
     &                    'DX', dx)


c$   Convert to km
        dx = dx / 1000.0

c$$$ Begin variable allocations

        allocate(lats(mix-1,mjx-1))
        allocate(lons(mix-1,mjx-1))
        allocate(fcstprobs(mix-1,mjx-1,mkx))
        allocate(obs(mix-1,mjx-1))
        allocate(prob_bins(nbins))
        allocate(fcst_freq(nbins))
        allocate(ob_hit_rate(nbins))

c$$$ Define probability bins

        prob_bins = (/ 10, 20, 30, 40, 50,
     &                  60, 70, 80, 90, 100 /)
        print*, prob_bins

c$$$ End variable allocations

         print*, mix-1,mjx-1,mkx,timo
         call get_variable3d(iunitprob,probvar,mix-1,
     &     mjx-1,mkx-1,timo,fcstprobs)
         call get_variable2d(iunitprob,"XLAT",mix-1,
     &     mjx-1,timo,lats)
         call get_variable2d(iunitprob,"XLONG",mix-1,
     &     mjx-1,timo,lons)
         call close_file(iunitprob)

c$$$ Open obs file and grab binary obs

        print*, "Observation variable", obvar
        call open_file(obfile,permissr,iunitob)
        rcode = nf_get_att_int(iunitob, nf_global,
     &                    'FHR', fcsthr)

        print*, "Fcst hr in ob file", fcsthr
        if (fcsthr .ne. fhr) then
          print*, "WARNING - OBSERVATION FILE FORECAST HOUR",
     &    " DOESN'T MATCH INPUT FORECAST HOUR!"
        endif
        print*, "Prob variable index", probvarind
c        rcode = nf_inq_varid(iunitob,obvar,obid)
c        rcode = nf_get_var_real(iunitob,obid,obs)
        call get_variable3d(iunitob,obvar,mix-1,
     &                          mjx-1,1,timo,obs)
        print*, "Obs Maximum", maxval(obs)
        call close_file(iunitob)

c$$$ Separate forecast probs into bins

        do p=1,nbins
         prob = prob_bins(p)
         do i=1,mix-1
          do j=1,mjx-1

c$$$ If forecast prob loc is in response box, separate into prob bin

            if ((lons(i,j) .ge. rbox_bounds(1)) .and.
     &          (lons(i,j) .lt. rbox_bounds(2)) .and.
     &          (lats(i,j) .ge. rbox_bounds(3)) .and.
     &          (lats(i,j) .lt. rbox_bounds(4))) then

c$$$ If forecast prob falls into bin, investigate further

            if ((abs(fcstprobs(i,j,probvarind)-prob) .le. 10)
     &          .and. (fcstprobs(i,j,probvarind) .lt. prob)) then
                   fcst_freq(p) = fcst_freq(p) + 1

c$$$ Handle edges

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
                       jbeg=j-15
                       jend=j+15
                    endif

c$$$ Check for points within neighborhood
                    obcnt = 0
                    do iii=ibeg,iend
                     do jjj=jbeg,jend
                       distz=(((iii-i)*dx)**2 + ((jjj-j)*dx)**2)**0.5
                       if (distz.le.nbrhd) then
                          if (obs(iii,jjj).gt.0) then
c                             print*, "Found prob that falls between..."
c                             print*, prob-10, "and", prob
c                             print*, "Prob in bin"
c                             print*, fcstprobs(i,j,probvarind)
c                             print*, "Prob loc", lats(i,j), lons(i,j)
c                             print*, "Distance of ob from prob loc"
c                             print*, distz
c                             print*, "Falls in neighborhood", nbrhd
c                             print*, "Ob", obs(iii,jjj)
                             obcnt = obcnt + 1
                          endif
                       endif
                     enddo
                    enddo
                    if (obcnt .gt. 0) then
                       ob_hit_rate(p) = ob_hit_rate(p) + 1
                    endif
           endif
           endif
          enddo
         enddo
         ob_hit_rate(p) = ob_hit_rate(p) / fcst_freq(p)
        enddo
        print*, "Max ob hit rate", MAXVAL(ob_hit_rate)
        print*, "Max fcst freq", MAXVAL(fcst_freq)

c$$$ Write results to file
        rcode = nf_create(outfile, nf_clobber, fidc)

        print*, "Created file", outfile

        bin = nbins
        print*, "Nbins", nbins

        rcode = nf_def_dim(fidc, 'nbins', nbins, bin)
        rcode = nf_def_var(fidc,'prob_bins',nf_float,
     &          1,bin,idprobbins)
        rcode = nf_put_att_text(fidc, idprobbins, 'description', 50,
     &                 'Probability Bins')
        rcode = nf_put_att_text(fidc, idprobbins,
     &                  'units', 20, 'percent')
        rcode = nf_put_att_real(fidc,idprobbins,'_FillValue',
     &                  nf_float,1,fmiss)

        rcode = nf_def_var(fidc,'fcst_frequency',nf_float,
     &          1,bin,idfcstfreq)
        rcode = nf_put_att_text(fidc, idfcstfreq, 'description', 50,
     &                 'Frequency of Probability Realizations')
        rcode = nf_put_att_text(fidc, idfcstfreq, 'units', 20,
     &                  'number of points')
        rcode = nf_put_att_real(fidc,idfcstfreq,'_FillValue',nf_float,1,
     &                     fmiss)

        rcode = nf_def_var(fidc,'ob_hit_rate',nf_float,
     &          1,bin,idobhitrate)
        rcode = nf_put_att_text(fidc, idobhitrate, 'description', 50,
     &                 'Probability Bins')
        rcode = nf_put_att_text(fidc, idobhitrate, 'units', 20,
     &                  'hits per total fcst')
        rcode = nf_put_att_real(fidc,idobhitrate,'_FillValue',
     &                  nf_float,1,fmiss)
        call close_file(fidc)

        call open_file(outfile, nf_write, fidc)

        call write_variable1d(fidc,'prob_bins',nbins,1,prob_bins)
        call write_variable1d(fidc,'fcst_frequency',nbins,1,fcst_freq)
        call write_variable1d(fidc,'ob_hit_rate',nbins,1,ob_hit_rate)
        call close_file(fidc)

        print*, "SUCCESSFUL COMPLETION OF RELIABILITY CALC"

        end
