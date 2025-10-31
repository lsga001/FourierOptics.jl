module Fourier

export fft2, ifft2

function fft2(g, delta)
  return fftshift(fft(fftshift(g)))*delta^2;
end

function ifft2(G, delta_f)
  N = size(G, 1);
  return ifftshift(ifft(ifftshift(G)))*(N*delta_f)^2;
end

function bluestein_scalar(
    Input::MonoBeam,
    z::Real,
    input_window::NumericalWindow,
    output_window::NumericalWindow)
  """Scalar diffraction computation method using Bluestein method. 
  Light: Science and Applications, DOI:10.1038/s41377-020-00362-z"""

  c = 299792458;

  nu = Input.frequency;
  lambda = c/nu;
  k = 2*pi/lambda;

  F0(x,y) = exp(im*k*z)/(im*lambda*z) * exp(im*k*(x^2 + y^2)/(2*z));
  F(u,v) = exp(im*k*(u^2 + v^2)/(2*z));

end
