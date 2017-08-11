function V = blk_tilyap(obj,x0)
% Linearizes the system about the state x0 and returns a candidate LyapunovFunction
% generated by solving the quadratic Lyapunov equation.
% @retval V an msspoly representation of the quadratic Lyapunov candidate 

if (~isTI(obj)) error('I don''t know that this system is time invariant.  Set the TI flags and rerun this method if you believe the system to be.'); end
if (~isCT(obj)) error('only handle CT case so far'); end

nX = getNumStates(obj);
typecheck(x0,'Point');
x0 = double(x0.inFrame(obj.getStateFrame)); 

u0 = zeros(getNumInputs(obj),1);

tol = 1e-10;
[A,B,C,D,xdot0] = linearize(obj,0,x0,u0);
if (any(abs(xdot0)>tol))
  xdot0
  error('f(x0) is not a fixed point');
end

S=.6228*eye(2);
frame=CoordinateFrame('LyapunovState',nX,'x');
frame.addTransform(AffineTransform(frame,obj.getStateFrame,eye(length(x0)),+x0));
obj.getStateFrame.addTransform(AffineTransform(obj.getStateFrame,frame,eye(length(x0)),-x0));

V = QuadraticLyapunovFunction(frame,S);

end