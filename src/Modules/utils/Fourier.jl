module Fourier

import FromFile: @from
@from "../Bare.jl" using Bare

export fft2, 
       ifft2, 
       chirpz, 
       czt, 
       fft_bluestein, 
       fft2_bluestein,
       bluestein_scalar

function fft2(g, delta)
  return fftshift(fft(fftshift(g)))*delta^2;
end

function ifft2(G, delta_f)
  N = size(G, 1);
  return ifftshift(ifft(ifftshift(G)))*(N*delta_f)^2;
end

function czt(x, A, W, M)
  # simple algorithm
  k = range(0, length=M);
  X = zeros(M);
  z = A*W^(-k);
  for n in range(0, length=N)
    X += x[n] * z^(-n);
  end
end

function chirpz(x,A,W,M)
  L = int(2^ceil(log2(M+N-1)));
  n = range(0,length=N);
  y = A^(-n) * W^(n^2/2) * x;
  Y = fft(y, L);

  v = zeros(L);
  v[1:M] = W^(-n[1:M]^2/2);
  v[L-N+1:end] = W^(-n[N-1:0:end]^2/2);
  V = fft(v);

  g = ifft(V*Y)[1:M];
  k = range(0, M);
  g *= W^(k^2/2);
  return g;
end

function fft_bluestein(x,f1,f2,fs,M;axis=0)
  """
  Based almost 1 to 1 on the code of Rafael Fuente in his Github 
  repository "diffractsim".
  """
  phi1 = 2*pi*f1/fs;
  phi2 = 2*pi*f2/fs;

  A0 = 1; W0 = 1;

  A = A0*exp(im * phi1);
  W = W0*exp(-im * (phi2-phi1)/(M-1)); # M-1 since we count segments, not points
  #return X = czt(x, A, W, M);
  return X = chirpz(x, A, W, M);
end

function fft2_bluestein(U, fx1, fx2, fxs, fy1, fy2, fys)
  (Nx, Ny) = size(U);
  return fft_bluestein(fft_bluestein(U, fx1, fx2, fxs, Nx, axis=0), fy1, fy2, fys, Ny, axis=1)
end

function bluestein_scalar(
    Input::MonoBeam,
    z::Real,
    input_window::NumericalWindow,
    output_window::NumericalWindow)
  """Scalar diffraction computation method using Bluestein method. 
  Light: Science and Applications, DOI:10.1038/s41377-020-00362-z"""
  # The basic algorithm can be described as E = F0 x fft2(E0 x F)

  c = 299792458;

  nu = Input.frequency;
  lambda = c/nu;
  k = 2*pi/lambda;

  xv = input_window.xv;
  yv = input_window.yv;
  uv = output_window.xv;
  vv = output_window.yv;

  F0(x,y) = exp(im*k*z)/(im*lambda*z) * exp(im*k*(x^2 + y^2)/(2*z));
  F(u,v) = exp(im*k*(u^2 + v^2)/(2*z));

  dx = xv[2] - xv[1];
  dy = yv[2] - yv[1];
  fxs = lambda*z/dx;
  fys = lambda*z/dy;

  fx1 = uv[1]+fxs/2;
  fx2 = uv[end]+fxs/2;
  fy1 = vv[1]+fxy/2;
  fy2 = vv[end]+fxy/2;

  return F0.(xv, yv') .* fft2_bluestein(Input .* F.(uv, vv'),fx1,fx2,fxs,fy1,fy2,fys);
end

end
