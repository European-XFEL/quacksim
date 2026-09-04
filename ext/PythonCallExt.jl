module PythonCallExt

import QuackSim
import PythonCall: Py, pyconvert

QuackSim.convert_Et(x::Py) = pyconvert(Float64, x)

end
