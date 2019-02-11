# things defined here would be defined using the pre-processor in C

global const ENABLE_ASSERTS = true

macro asserts_enabled(expr1)

  if ENABLE_ASSERTS
    return quote
      $(esc(expr1))
    end
  else
    return nothing
  end

end  # end macro

global const HAVE_OPTIM = false # whether Optim is installed or not

"""
  Function throws an error with a helpful message if Optim is required but
  not available (based on the `HAVE_OPTIM` constant).
"""
function requireOptim()
  if !HAVE_OPTIM
    error("SummationByParts is not configured to use Optim, cannot build minimum condition number operator")
  end

  return nothing
end


