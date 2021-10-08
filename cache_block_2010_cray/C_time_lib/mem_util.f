       subroutine print_hwm
       parameter (lun=99)
       character*80 line
       integer(kind=8) hwm,rss
       open(lun,file='/proc/self/status')
       hwm=0
       do while(.true.)
         read(lun,'(a)',end=99) line
         if (line(1:6).eq.'VmHWM:') read(line(8:80),*) hwm
         if (line(1:6).eq.'VmRSS:') read(line(8:80),*) rss
       enddo
   99  close(lun)
       print *,'HWM = ',hwm,'kB =',hwm*1024,'bytes'
       print *,'RSS = ',rss,'kB =',rss*1024,'bytes'
       return
       end

