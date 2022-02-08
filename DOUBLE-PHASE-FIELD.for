*     *************************************************************************************************
*
*     Common variable module, which defines some common variables 
*     and allocates memory
*
*     *************************************************************************************************
*      
      module vars_module
          parameter (NumEle= 100000)
          real*8,save :: allD(NumEle),allH(NumEle)
          real*8,save :: allft(NumEle)
      end module

*     *************************************************************************************************
*
*     VUMAT for Visulization 
*
*     *************************************************************************************************     
*      
      subroutine vumat(
* Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
* Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
*
      use vars_module
      include 'vaba_param.inc'
*
* All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),
     4  stretchOld(*), defgradOld(*),
     5  fieldOld(*), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(*), fieldNew(*),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
*
      character*80 cmname
*
      parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     2  threeHalfs = 1.5 )
      
      dimension e_voigt(6),e_tensor(3,3),e_prin(3),e_dir(3,3)
      dimension alphai(3)
      dimension s_tensor(3,3), s_prin(3,3)

      integer jElemUid
      
      real*8  lc,lb, Gf, ea, ft 
      real*8  p, a1, a2, a3
      real*8  phase, dalpha, dphase
      real*8  omega,domega
      real*8  hist, hist_threshold
      real*8  c0
*
*     ================================================================================================     *
*     Statev(1)    Element Identification Number
*
*     ndir = 3
*     nshr = 1
*     ================================================================================================     *
*     Parameter initilization
      xnu   = 0.0
*      
      ea   =  props(1)      ! props(1) -- Young's modulus
      ft   =  props(2)       ! props(3) -- failure strength
      Gf   =  props(3)      ! props(4) -- fracture energy
      lb   =  props(4)      ! props(5) -- length scale  
      thk  =  1.0d0       ! props(6) -- thickness
*
      c0   =  3.1415926535897932384626433832d0
*
      twomu  = ea / ( one + xnu )
      alamda = twomu * xnu / ( one - two * xnu )
*     ================================================================================================     *
      if ( stepTime .eq. zero ) then
          do k = 1, nblock
              trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
*
              stressNew(k,1) = stressOld(k,1) 
     *                         + twomu * strainInc(k,1) + alamda * trace
              stressNew(k,2) = stressOld(k,2) 
     *                         + twomu * strainInc(k,2) + alamda * trace
              stressNew(k,3) = stressOld(k,3) 
     *                         + twomu * strainInc(k,3) + alamda * trace
              stressNew(k,4) = stressOld(k,4) + twomu * strainInc(k,4)
          end do
      end if
*     ================================================================================================     *
      if ( stepTime .gt. zero ) then
*
      do  k = 1, nblock
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
          jElemUid      =  stateOld(k,8)
*        
          a1   =  4.d0/(c0*lb)*ea*Gf/(ft*ft)
*         concrete softening
          p  =  2.0d0
          a2 =  1.3868d0
          a3 =  0.6567d0
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         Elastic strain
          stateNew(k,1) = stateOld(k,1) + strainInc(k,1)
          stateNew(k,2) = stateOld(k,2) + strainInc(k,2)
          stateNew(k,3) = stateOld(k,3) + strainInc(k,3)
          stateNew(k,4) = stateOld(k,4) + strainInc(k,4)
          stateNew(k,5) = 0.0d0
          stateNew(k,6) = 0.0d0
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
          e_voigt = stateNew(k,1:6)

*         The transformation of elastic strain from Voigt form to tensor form
          call voigt_convection(e_voigt, e_tensor,.false.,.true.)
          
*         Eigenvalue decomposition of strain
          call eigen3D(e_tensor, e_prin, e_dir)
*
          e_tr   = e_prin(1) + e_prin(2) + e_prin(3)
          e1plus = max(e_prin(1),0.0)
          e2plus = max(e_prin(2),0.0)
          e3plus = max(e_prin(3),0.0)      

          stressPower = half*(alamda*e_tr**two 
     1           + twomu*(e1plus**two + e2plus**two + e3plus**two) )
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *

          phase = allD(jElemUid) 
          
          if (props(1) .eq. 2) then
              phase = 0.0
          end if
*          
          call energeticFunc(omega,domega,phase,a1,a2,a3,p)
*
          do i=1,3
              alphai(i) = zero
              if (e_prin(i) .gt. zero) alphai(i)=one
          enddo
*          
          alpha = one
*

          
          if (e_prin(1) .gt. 0.0d0) then
              sigma1 =  twomu * e_prin(1) * omega
     1            + alamda * e_tr * omega
          else 
               sigma1 =  twomu * e_prin(1) 
     1            + alamda * e_tr 
          end if
*          
          if (e_prin(2) .gt. 0.0d0) then
              sigma2 =  twomu * e_prin(2) * omega
     1            + alamda * e_tr * omega
          else 
               sigma2 =  twomu * e_prin(2) 
     1            + alamda * e_tr 
          end if
          

          
          s_prin = 0.d0
          s_prin(1,1) = sigma1
          s_prin(2,2) = sigma2
          s_prin(3,3) = sigma3

          s_tensor = matmul(matmul(e_dir,s_prin),transpose(e_dir))
          
          call voigt_convection(stressNew(k,:), s_tensor,.true.,.true.)
          
          stressNew(k,1) = stressNew(k,1)
          stressNew(k,2) = stressNew(k,2)
          stressNew(k,3) = 0.0d0
          stressNew(k,4) = stressNew(k,4)
          stressNew(k,5) = 0.0d0
          stressNew(k,6) = 0.0d0
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         stress power      
          stateNew(k,11) =  stressPower
*         histroy variable, H  
          stateNew(k,14) =  max(stateOld(k,14), stateNew(k,11))
*           
          allH(jElemUid) = stateNew(k,14) 
*          
          stateNew(k,15) =  phase
          
          stateNew(k,9) =   sigma1
          stateNew(k,10) =  sigma2
          stateNew(k,12) =  omega
          stateNew(k,13) =  allft(jElemUid)
*     ================================================================================================     *
      end do
      end if
*
      return
      end

      
      
*     *************************************************************************************************
*
*     VUEL for the phase field  
*     Activated DOFs: 11
*
*     *************************************************************************************************
      subroutine vuel(
     *     nblock,
*          to be defined
     *     rhs,amass,dtimeStable,
     *     svars,nsvars,
     *     energy,
*          
     *     nnode,ndofel,
     *     props,nprops,
     *     jprops,njprops,
     *     coords,ncrd,
     *     u,du,v,a,
     *     jtype,jelem,
     *     time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,
     *     dMassScaleFactor,
     *     predef,npredef,
     *     jdltyp,adlmag)
* 
      use vars_module
      include 'vaba_param.inc'

      parameter ( zero = 0.d0, half = 0.5d0, one = 1.d0, two=2.d0 )

*     operation code
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)

*     flags
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)

*     time
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)

*     procedure flags
      parameter ( jDynExplicit = 17 )

*     energies 
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)

*     predefined variables
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)    

*     indexing in a 3-long vector

      parameter (factorStable = 0.99d0)
      
      dimension rhs(nblock,ndofel), amass(nblock,ndofel,ndofel),
     *     dtimeStable(nblock),
     *     svars(nblock,nsvars), energy(nblock,nElEnergy),
     *     props(nprops), jprops(njprops),
     *     jelem(nblock), time(nTime), lflags(nFlags),
     *     coords(nblock,nnode,ncrd), u(nblock,ndofel),
     *     du(nblock,ndofel), v(nblock,ndofel), a(nblock, ndofel),
     *     dMassScaleFactor(nblock),
     *     predef(nblock, nnode, npredef, nPred), adlmag(nblock)
*
      dimension UD(4),DUD(4)
*
      integer i,j,k
      integer ngp
*      
      real*8  lc,lb, Gf, ea, ft 
      real*8  p, a1, a2, a3
      real*8  phase, dalpha, dphase
      real*8  omega,domega
      real*8  hist, hist_threshold
      real*8  c0
      
      parameter (ngp = 2)
      
      dimension gp(ngp), gw(ngp)
      dimension bd(2,4) 
      dimension coords_ele(2, 4)
      dimension shN(4), dn_xieta(2,4), dn_x(4), dn_y(4)
      dimension rd(4)

*     ================================================================================================     *
*     parameter initilization
*
      eta = 1.0d-5
*
      ea   =  props(1)      ! props(1) -- Young's modulus
      ft   =  props(2)      ! props(3) -- failure strength
      Gf   =  props(3)      ! props(4) -- fracture energy
      lb   =  props(4)      ! props(5) -- length scale 
      
      thk  =  1.0d0       ! props(6) -- thickness
*
      c0   =  3.1415926535897932384626433832d0
*
      WF = Gf/2.0/lb      
      WO = zero  
*
      gp = (/ -1.d0, 1.d0 /) / dsqrt(3.d0)
      gw = (/  1.d0, 1.d0 /)
*
      if ( lflags(iOpCode).eq.jMassCalc ) then
*     ================================================================================================     *
          do kblock = 1, nblock
              amass(kblock,1,1) = eta *10.0d0
              amass(kblock,2,2) = eta *10.0d0
              amass(kblock,3,3) = eta *10.0d0
              amass(kblock,4,4) = eta *10.0d0            
          end do
*     ================================================================================================     *
      else if ( lflags(iOpCode) .eq. jIntForceAndDtStable) then
*     ================================================================================================     *
          do kblock = 1, nblock
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *          
              a1   =  4.d0/(c0*lb)*ea*Gf/(ft*ft)
*             concrete softening
              p  =  2.0d0
              a2 =  1.3868d0
              a3 =  0.6567d0 
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *          
*             coords_ele(ncrd, nnode) = coords(kblock,nnode,ncrd)
*             nnode = 4; ncrd = 3
              coords_ele(1, 1) = coords(kblock,1,1)
              coords_ele(2, 1) = coords(kblock,1,2)
              coords_ele(1, 2) = coords(kblock,2,1)
              coords_ele(2, 2) = coords(kblock,2,2)
              coords_ele(1, 3) = coords(kblock,3,1)
              coords_ele(2, 3) = coords(kblock,3,2)
              coords_ele(1, 4) = coords(kblock,4,1)
              coords_ele(2, 4) = coords(kblock,4,2)
              
              do i = 1, nnode
                  UD(i)  = U(kblock,i)
                  DUD(i) = DU(kblock,i) 
              end do
              
              hist = allH(jelem(kblock) - NumEle)
              if(hist .lt. zero)  hist = zero 
              
              call geometricFunc(dalpha,0.0d0) 
              call energeticFunc(omega,domega,0.0d0,a1,a2,a3,p)  
              hist_0 = - Gf * dalpha / domega / c0 / lb
              
              hist = max(hist, hist_0)
*              
              rd = 0.0d0
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *          
              do i = 1, ngp
                do j = 1, ngp
*                    
                  call b_matrix(shN,bd,det_jacb,coords_ele,gp(i),gp(j))
*                      
                  phase  = dot_product(shN,UD)                  
                  
                  call geometricFunc(dalpha,phase)                                ! geometric function
                  call energeticFunc(omega,domega,phase,a1,a2,a3,p)               ! energetic function
*
                  phi_source  = domega *  hist + Gf/(c0*lb)*dalpha
                  
                  dvol=  gw(i)*gw(j)*det_jacb*thk
*
                  rd  =  rd  - dvol*(phi_source*shN + 2.d0*lb*Gf/c0
     &                    *  matmul(transpose(bd), matmul(bd, UD)))    
                  
                  rd  =  rd  - dvol*(phi_source*shN )    
*
                  end do
                end do
*
              rhs(kblock,1) =  - rd(1)
              rhs(kblock,2) =  - rd(2)
              rhs(kblock,3) =  - rd(3)
              rhs(kblock,4) =  - rd(4)
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *          
              if (time(2) .le. 2.0) then 
*              if (hist .le. hist_0) then 
                  rhs(kblock,1) = 0.0
                  rhs(kblock,2) = 0.0
                  rhs(kblock,3) = 0.0
                  rhs(kblock,4) = 0.0
              end if
              
              allD(jelem(kblock) - NumEle)  = phase    
              allft(jelem(kblock)- NumEle)  = ft
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *          
          end do
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
      end if
*     ================================================================================================     *
      return
      end      
      
      
      
*     *************************************************************************************************
*
*     VUSDFLD for converting the state variables in VUMAT into global variables  
*     
*     *************************************************************************************************
*       
      subroutine vusdfld(
* Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElemUid, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
* Write only -
     *   stateNew, field )
*
      include 'vaba_param.inc'
*
      dimension props(nprops),
     *          jElemUid(nblock), coordMp(nblock, *), 
     *          direct(nblock, 3, 3), T(nblock,3,3), 
     *          charLength(nblock),
     *          stateOld(nblock, nstatev), 
     *          stateNew(nblock, nstatev),
     *          field(nblock, nfieldv)
      character*80 cmname
*
      character*3 cData(maxblk)
      dimension jData(maxblk)
      dimension eqps(maxblk)
*
      parameter ( zero = 0.d0 )
*     ================================================================================================     *
      do k = 1, nblock
            stateOld(k,8)   =  jElemUid(k)
            stateNew(k,8)   =  jElemUid(k)
      end do
*     ================================================================================================     *
      return
      end subroutine
     
      
C     *************************************************************************************************
C
C     energetic degradation function  omega  
C     
C     *************************************************************************************************
C        
      subroutine energeticFunc(omega,domega,phi,a1,a2,a3,p)
*
      include 'vaba_param.inc'
      
      real*8:: omega, domega, phi
      real*8:: fac1, dfac1, ddfac1, fac2, dfac2, ddfac2
      real*8:: p, a1, a2, a3
*     ================================================================================================     *
*      
      fac1    =  (1.d0 - phi)**p
      dfac1   = - p*(1.d0 - phi)**(p - 1.d0); 
      ddfac1  =  p*(p - 1.d0)*(1.d0 - phi)**(p - 2.d0)
*        
      fac2   =  fac1   + a1*phi + a1*a2*phi**2.d0 + a1*a2*a3*phi**3.d0
      dfac2  =  dfac1  + a1 + 2.d0*a1*a2*phi + 3.d0*a1*a2*a3*phi**2.d0
      ddfac2  =  ddfac1 + 2.d0*a1*a2 + 6.d0*a1*a2*a3*phi
*        
      omega   =  fac1/fac2        
      domega  =  (dfac1*fac2  - fac1*dfac2)/(fac2**2.d0)
*       
*     ================================================================================================     *
      return
      end subroutine energeticFunc
      
!**********************************************************************************************************
!
      subroutine geometricFunc(dalpha,phase)
!
!**********************************************************************************************************
      include 'vaba_param.inc'
      real*8  dalpha, phase
*     ================================================================================================     *
*        
      dalpha  = 2.d0 - 2.d0*phase
*     ================================================================================================     *
        
      return 
      end subroutine geometricFunc  

      
*     *************************************************************************************************
*
*     shape function 
*     
*     *************************************************************************************************
*        
      subroutine shapefuc(shN, dn_xieta, xi, eta)
          
      include 'vaba_param.inc'
      dimension shN(4), dn_xieta(2, 4) 
*     ================================================================================================     *
      shN(1) = 0.25d0*(1.d0 - xi)*(1.d0 - eta)
      shN(2) = 0.25d0*(1.d0 + xi)*(1.d0 - eta)
      shN(3) = 0.25d0*(1.d0 + xi)*(1.d0 + eta)
      shN(4) = 0.25d0*(1.d0 - xi)*(1.d0 + eta)
*          
      dn_xieta(1, 1) = -0.25d0*(1.d0 - eta)
      dn_xieta(1, 2) =  0.25d0*(1.d0 - eta)
      dn_xieta(1, 3) =  0.25d0*(1.d0 + eta)
      dn_xieta(1, 4) = -0.25d0*(1.d0 + eta)
*            
      dn_xieta(2, 1) = -0.25d0*(1.d0 - xi)
      dn_xieta(2, 2) = -0.25d0*(1.d0 + xi)
      dn_xieta(2, 3) =  0.25d0*(1.d0 + xi)
      dn_xieta(2, 4) =  0.25d0*(1.d0 - xi)
*     ================================================================================================     *
      return 
      end subroutine shapefuc 
      

*     *************************************************************************************************
*
*     B matrix
*     
*     *************************************************************************************************
*        
      subroutine b_matrix(shN,bd,det_jacb,coords,xi,eta)
      
      include 'vaba_param.inc'
      dimension bd(2,4) 
      dimension djacb(2,2), dinv_jacb(2,2), coords(2, 4)
      dimension shN(4), dn_xieta(2,4), dn_x(4), dn_y(4)
*     ================================================================================================     *
*     shape functions 
      call shapefuc(shN, dn_xieta, xi, eta)
            
*     jacob matrix
      djacb = matmul(dn_xieta, transpose(coords))            
      det_jacb       = djacb(1,1)*djacb(2,2) - djacb(1,2)*djacb(2,1)
      dinv_jacb(1, 1) =   djacb(2, 2)
      dinv_jacb(1, 2) = - djacb(1, 2)
      dinv_jacb(2, 1) = - djacb(2, 1)
      dinv_jacb(2, 2) =   djacb(1, 1)
      dinv_jacb       =   1.d0 / det_jacb * dinv_jacb    
        
*     initialize varibles
      do i = 1,4
         dn_x(i) = dinv_jacb(1,1)*dn_xieta(1,i) + 
     &             dinv_jacb(1,2)*dn_xieta(2,i)
         dn_y(i) = dinv_jacb(2,1)*dn_xieta(1,i) + 
     &             dinv_jacb(2,2)*dn_xieta(2,i)
      end do
           
*     B matrix for damage
      do j = 1,4
         bd(1,j) = dn_x(j)
         bd(2,j) = dn_y(j)
      end do
*     ================================================================================================     *
      return
      end subroutine b_matrix

      
      
c ----------------------------------------------------------------------
C
c Convert between voigt notation and tensor form. 
C it only works for 3D condition. 
C
c ----------------------------------------------------------------------
      subroutine voigt_convection(A_voigt, A, convert2voigt, kinetic)
      ! Checked.
      integer i,j
      real, dimension(3,3):: A
      real A_voigt(6)
      logical convert2voigt, kinetic

      if (convert2voigt) then 
          ! Convect to the voigt notation.
          A_voigt(:) = 0
          do i = 1, 3
              A_voigt(i) = A(i,i)        
          end do

          if (kinetic) then
              A_voigt(4) = A(1,2)
              A_voigt(5) = A(2,3)
              A_voigt(6) = A(1,3)
          else  
              A_voigt(4) = 2.0*A(1,2)
              A_voigt(5) = 2.0*A(1,3)
              A_voigt(6) = 2.0*A(2,3)
          end if
      else 
          ! Convect to the tensor form.
          A(:,:) = 0
          do i = 1, 3
                A(i,i) = A_voigt(i)        
          end do
          if (kinetic) then
                A(1,2) =  A_voigt(4) 
                A(2,3) =  A_voigt(5) 
                A(1,3) =  A_voigt(6) 
          else 
                A(1,2) = 0.5 * A_voigt(4) 
                A(1,3) = 0.5 * A_voigt(5) 
                A(2,3) = 0.5 * A_voigt(6) 
          end if 
          do i = 1,2
              do j = i+1,3
                  A(j,i) = A(i,j)
              end do
          end do
      end if
      
      return
      end subroutine
      
      
      
c ----------------------------------------------------------------------
C
C    eigen3D = determine eigenvalues/vectors for 3x3 symmetric matrix
C    
C    arguments description
C    ---------------------
C    input:
C      u(3,3) : matrix with initial values (only upper half used)
C    
C    outputs:
C      d(3) : eigenvalues associated with columns of v
C      v(3,3) : matrix of the eigenvectors (by column)
C
c ----------------------------------------------------------------------

      subroutine eigen3D(u, d, v)

      real, dimension(3,3), intent(in) :: u
      real, dimension(3), intent(out) :: d
      real, dimension(3,3) :: v
      ! ====================================
      ! local variables
      ! ==============
      integer :: rot, its
      real :: g, h, aij, sm, thresh, t, c, s, tau

      real :: a(3), b(3), z(3)

      ! loop index variable
      integer :: i, j, k
      
      parameter ( zero = 0.d0, one = 1.d0 , two = 2.d0)
      ! ====================================

      ! initialize
      d(:)= zero
      v(:,:)= zero

      ! copy value of u to v
      v(1:3,1:3)= u(1:3,1:3)

      ! transfer arrays into 1D arrays
      a(1) = v(1,2)
      a(2) = v(2,3)
      a(3) = v(1,3)

      do i = 1, 3
         d(i) = v(i,i)
         b(i) = d(i)
         z(i) = zero

         do j = 1,3
            v(i,j) = zero
         end do ! j

         v(i,i) = one
      end do ! i

      ! check for diagonal case
      sm = abs(a(1)) + abs(a(2)) + abs(a(3))
      g  = abs(d(1)) + abs(d(2)) + abs(d(3))

      if (sm .lt. 1.0e-13*g) return

      rot = 0
      do its = 1, 50

         ! set convergence test and threshold
         sm = abs(a(1)) + abs(a(2)) + abs(a(3))
         if ( sm == zero ) return

         if( its < 4 ) then
            thresh = 0.011d0*sm
         else
            thresh = zero
         end if

         ! perform sweeps for rotation
         do i = 1, 3
            j = mod(i,3) + 1
            k = mod(j,3) + 1

            aij = a(i)
            g = 100.d0 * abs(aij)

          if((abs(d(i))+g/=abs(d(i))).or.(abs(d(j))+g/=abs(d(j)))) then

               if( abs(aij) > thresh ) then

                  a(i) = zero
                  h = d(j) - d(i)

                  if( abs(h)+g == abs(h) ) then
                     t=aij / h
                  else
                     t=sign(two,h/aij)/(abs(h/aij)+sqrt(4.0+(h/aij)**2))
                  end if

                  ! set rotation parameters
                  c = one/sqrt(one+t*t)
                  s = t * c
                  tau = s / (one+c)

                  ! rotate diagonal terms
                  h = t * aij
                  z(i) = z(i) - h
                  z(j) = z(j) + h
                  d(i) = d(i) - h
                  d(j) = d(j) + h

                  ! rotate off-diagonal terms
                  h = a(j)
                  g = a(k)
                  a(j) = h + s*(g - h*tau)
                  a(k) = g - s*(h + g*tau)

                  ! rotate eigenvectors
                  do k = 1, 3
                     g = v(k,i)
                     h = v(k,j)
                     v(k,i) = g - s*(h + g*tau)
                     v(k,j) = h + s*(g - h*tau)
                  end do ! k

                  rot = rot + 1

               end if

            else

               a(i) = zero

            end if

         end do ! i

         ! update diagonal terms
         do i = 1, 3
            b(i) = b(i) + z(i)
            d(i) = b(i)
            z(i) = zero
         end do ! i

      end do ! its

      return
      end subroutine