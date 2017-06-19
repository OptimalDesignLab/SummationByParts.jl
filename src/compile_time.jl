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
