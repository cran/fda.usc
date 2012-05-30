      Subroutine aijr(n,inprod,Aijr_vec)
      
	  ! Arguments
      Double Precision inprod(n,n), Aijr_vec((n*n*n+n*n)/2)
      Integer n
      
	  ! Local variables
      Double Precision num, den, coc, aux1, aux2
      Integer i,j,r
      Real, parameter :: Pi = 3.1415926536
	   
	  
	  ! Calculus of Aijr
      do r=1,n
      do i=1,n
      do j=1,i
	  
        if(i==j) then
          if(i==r) then 
            Aijr_vec((r-1)*(n*n+n)/2+(i*(i-1))/2+j)=2*PI
          else 
            Aijr_vec((r-1)*(n*n+n)/2+(i*(i-1))/2+j)=PI
          end if
        else
          if((i==r) .or. (j==r)) then
            Aijr_vec((r-1)*(n*n+n)/2+(i*(i-1))/2+j)=PI
          else
            num=inprod(i,j)-inprod(i,r)-inprod(r,j)+inprod(r,r)
            aux1=sqrt(inprod(i,i)-2*inprod(i,r)+inprod(r,r))
            aux2=sqrt(inprod(j,j)-2*inprod(j,r)+inprod(r,r))
            den=aux1*aux2
            coc=num/den
            if(coc<-1) then
              coc=-1
            else if(coc>1) then
              coc=1
            end if
            Aijr_vec((r-1)*(n*n+n)/2+(i*(i-1))/2+j)=abs(PI-acos(coc))
          end if
          !Aijr_vec(j,i,r)=Aijr_vec(i,j,r)
        end if
	  
      end do 			
      end do
      end do
      	   
      Return
      End