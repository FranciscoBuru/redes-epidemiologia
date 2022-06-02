#Pkg.add(PackageSpec(name="PyCall", rev="master"))
#Pkg.build("PyCall")
using PyCall

function hit_and_run(puntos::Integer = 2)
  py"""
  from hitAndRun import hit_and_run
  import numpy as np
  def puntos(num = 10):
    return np.nan_to_num(hit_and_run(num), nan = 1.0)
  """
  return py"puntos"(puntos)
end

print(hit_and_run(45))
