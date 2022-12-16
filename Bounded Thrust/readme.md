Bounded Thrust Accleration to 5*1e-4 m/s^2        
### 1. BT_step1              
  
find initial costate by steepest gradient method (1st order)                   
#### trouble:           
Unbounded Thrust 문제에서 찾은 costate (initial guess) 에서 처음 한 번 gradient 반대 방향으로 costate update 시켜주면 error가 커짐        
그 이후 error가 줄긴 하지만 차이가 미미함        
이상한 local minima에 빠진 듯

 
### 2. UBT               
Unbounded Thrust 문제에서 찾은 costate        
BT_step1의 costate inital guess로 쓰인다.              

### 3. BT_eqn         
Bounded Thrust 문제의 state & costate adjoint vector 적분 
