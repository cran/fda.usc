      Subroutine pcvm_statistic(n,Aijr_vec,residuals,statistic)
      
	  ! Arguments
      Double Precision statistic, Aijr_vec((n*n*n+n*n)/2), residuals(n)
      Integer n
      
	  ! Local variables
      Double Precision statistic1, statistic2, aux1, aux2
      Integer i,j,r   
	  
      statistic1=0
      statistic2=0
	  
	  ! Calculus of Aijr
      do r=1,n
      do i=1,n
      do j=1,i
       
       aux1=Aijr_vec((r-1)*(n*n+n)/2+(i*(i-1))/2+j)
       aux2=residuals(i)*residuals(j)
	   
       if(i==j) then
         statistic1=statistic1+aux1*aux2
       else
         statistic2=statistic2+aux1*aux2
       end if
	   	  
      end do 			
      end do
      end do
      	 
      statistic=statistic1+2*statistic2		 
	  
      Return
      End