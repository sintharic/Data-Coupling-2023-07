rey() {
  rm -f currentCalc.x
  cp ~/01src/currentCalc.x . 
  ./currentCalc.x "$@" > output.log &
}