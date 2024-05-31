assignInNamespace(
  "cli_blue_bullet",
  function (msg, print = TRUE) 
  {
      symbol <- targets:::cli_symbol_bullet_blue
      msg <- cli::col_blue(paste(symbol, msg))
      targets:::if_any(print, message(msg), msg)
  },
  ns = "targets"
)

if (F)
assignInNamespace(
  "cli_green_check",
function (msg, print = TRUE) 
{
    symbol <- targets:::cli_symbol_tick_green
    msg <- cli::col_blue(paste(symbol, msg))
    targets:::if_any(print, message(msg), msg)
},
ns = "targets"
)