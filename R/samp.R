samp <- function(knots, max) {
  if((knots == 0) | (knots == 1)) {return(1)} # having 0 or 1 knots must auto defer to birth
  if(knots == max) {sample(2:3, 1)} #at max_knot knot capacity, can only delete or change
  sample(3, 1)
  # 1 = BIRTH
  # 2 = DEATH
  # 3 = CHANGE
}
