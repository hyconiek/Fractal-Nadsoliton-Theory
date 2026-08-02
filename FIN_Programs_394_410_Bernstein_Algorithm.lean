import Std

/- P394: exact reflection of the complete depth-14 dyadic Bernstein
   subdivision algorithm.  Every cell uses a common positive denominator.
   One degree-11 split multiplies that denominator by 2^11. -/
namespace FINPrograms394To410Bernstein

def degree : Nat := 11
def depth : Nat := 14
def initialDenominator : Int := (28544954956131505450700927817744492475466705001513737330793735 : Int)
def initialNumerators : List Int := [
  (-17608035311570273326721990545666495915798036573079500162097092 : Int),
  (-135603318270838358135208498520498224151290362456984784759971780 : Int),
  (637959137993784129012225963874383825458128573134989497722183740 : Int),
  (-1702154343499180893828421226443100499882529614751891558442229700 : Int),
  (3180994873242910514425865943013909753574550136689975166819082300 : Int),
  (-4354651047557224726775218360114431683814134452211543645771816900 : Int),
  (4318374090595628182208856005668869556868756497835702837841469500 : Int),
  (-2956518923803658677909561949963072491423699817797656781517250500 : Int),
  (1212824592942340589249930657960610372245366988363965927001661500 : Int),
  (-183235308335264173780008218358351775213759525265756408099921860 : Int),
  (-75624098027535634180833892169169846864726475519775166977040324 : Int),
  (-20360944239790424575424925125936933551647811432898500 : Int)
]

def choose : Nat -> Nat -> Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

def leftCoefficient (cell : List Int) (j : Nat) : Int :=
  ((List.range (j + 1)).foldl
    (fun total i => total + (choose j i : Int) * cell.getD i 0)
    0) * (2 : Int) ^ (degree - j)

def rightCoefficient (cell : List Int) (j : Nat) : Int :=
  ((List.range (degree - j + 1)).foldl
    (fun total i => total + (choose (degree - j) i : Int) *
      cell.getD (j + i) 0)
    0) * (2 : Int) ^ j

def splitCell (cell : List Int) : List (List Int) :=
  [(List.range (degree + 1)).map (leftCoefficient cell),
   (List.range (degree + 1)).map (rightCoefficient cell)]

def cellSafe (denominator : Int) (cell : List Int) : Bool :=
  cell.all (fun value => decide (-denominator <= value && value <= 0))

def certify : Nat -> Int -> List Int -> Bool
  | 0, denominator, cell => cellSafe denominator cell
  | remaining + 1, denominator, cell =>
      if cellSafe denominator cell then true
      else
        let children := splitCell cell
        let nextDenominator := denominator * (2 : Int) ^ degree
        certify remaining nextDenominator (children.getD 0 []) &&
          certify remaining nextDenominator (children.getD 1 [])

def nodesVisited : Nat -> Int -> List Int -> Nat
  | 0, _, _ => 1
  | remaining + 1, denominator, cell =>
      if cellSafe denominator cell then 1
      else
        let children := splitCell cell
        let nextDenominator := denominator * (2 : Int) ^ degree
        1 + nodesVisited remaining nextDenominator (children.getD 0 []) +
          nodesVisited remaining nextDenominator (children.getD 1 [])

def certified : Bool := certify depth initialDenominator initialNumerators

theorem complete_subdivision_certificate : certified = true := by
  native_decide

theorem reflected_adaptive_node_count :
    nodesVisited depth initialDenominator initialNumerators = 147 := by
  native_decide

end FINPrograms394To410Bernstein
