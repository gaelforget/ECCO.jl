module ECCOEnzymeExt

import ECCO: _autodiff_Reverse, _Active, _Duplicated
import Enzyme: autodiff, Reverse, Active, Duplicated

_Duplicated(args...; kwargs...) = Duplicated(args...; kwargs...)
_Active(args...; kwargs...) = Active(args...; kwargs...)
_autodiff_Reverse(args...; kwargs...) = autodiff(Reverse, args...; kwargs...)

end


