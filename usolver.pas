{$mode ObjFPC}{$H+}
unit uSolver;
interface
type TVector = array of real;
     {тип для передачи правых частей уравнения
      inv - аргумент функции
      outv - значение функции
      ExtraData - дополнительные данные }
     TRPFunc = procedure (const inv: TVector; var outv: TVector;
                          ExtraData: pointer = nil);

     {функция печати}
     TPrintProc = procedure (const inv: TVector; ExtraData: pointer = nil);


     {Решение системы ОДУ методом Рунге-Кутты 4-го порядка с фиксированным шагом
      ВХОД
      y0 - начальные условия
      rp - функция правых частей
      Step - шаг интегрирования
      exitcond - функция выхода
      nconstr - количество функция выхода
      TOL - точность выхода
      pf - функция печати
      ExtraData - указатель на дополнительные данные задачи
      РЕЗУЛЬТАТ
      Функция возвращает номер сработавшей функции выхода.
      y0 будет содержать значение фазовых координат в последней точке.
      }
     function RK4Fixed(var y0: TVector; rp: TRPFunc; Step: real;
                       exitcond: TRPFunc; nconstr: Integer; TOL: real = 1E-6;
                       pf: TPrintProc = nil; ExtraData: pointer = nil):Integer;

implementation

function RK4Fixed(var y0: TVector; rp: TRPFunc; Step: real; exitcond: TRPFunc;
  nconstr: Integer; TOL: real; pf: TPrintProc; ExtraData: pointer): Integer;
var tmp,         // временный вектор
    y1, y2: TVector; // следующая точка
    k: array of TVector;
    neqn: Integer; // число уравнений

procedure DoStep(y0, y1: TVector);
var i: Integer;
begin
 // k1
 rp(y0, k[0], ExtraData);

 //k2
 tmp[0]:=y0[0] + Step/2;
 for i:=1 to neqn do
  tmp[i]:=y0[i] + 0.5*Step*k[0, i];
 rp(tmp, k[1], ExtraData);

 //k3
 tmp[0]:=y0[0] + Step/2;
 for i:=1 to neqn do
  tmp[i]:=y0[i] + 0.5*Step*k[1, i];
 rp(tmp, k[2], ExtraData);

 //k4
 tmp[0]:=y0[0] + Step;
 for i:=1 to neqn do
  tmp[i]:=y0[i] + Step*k[2, i];
 rp(tmp, k[3], ExtraData);

 //y_k+1
 y1[0]:=y0[0] + Step;
 for i:=1 to neqn do
  y1[i]:=y0[i] + Step/6*(k[0, i] + 2*k[1, i] + 2*k[2, i] + k[3, i]);
end;

function CheckCond(const x1, x2: TVector): Integer;
var i: Integer;
begin
 Result:=-1;
 for i:=0 to High(x1) do
  if x1[i]*x2[i] <= 0 then begin
   Result:=i;
   Exit;
  end;
end;

var i: Integer;
    ext1,          // значение функций выхода в точке y0
    ext2,
    ext3: TVector; // значение функций выхода в точке y1
begin
 Result:=-1;

 // создание временных переменных
 neqn:=Length(y0) - 1;
 SetLength(tmp, neqn + 1);
 SetLength(y1, neqn + 1);
 SetLength(y2, neqn + 1);
 SetLength(k, 4, neqn + 1);
 SetLength(ext1, nconstr);
 SetLength(ext2, nconstr);
 SetLength(ext3, nconstr);

 exitcond(y0, ext1, ExtraData);
 if Assigned(pf) then
  pf(y0, ExtraData);

 while True do begin
  DoStep(y0, y1);

  exitcond(y1, ext2, ExtraData);
  Result:=CheckCond(ext1, ext2);
  if Result >= 0 then
   break;

  if Assigned(pf) then
   pf(y1, ExtraData);

  y0:=Copy(y1);
  ext1:=Copy(ext2);
 end;

 // уточнение решения методом дихотомии
 while (abs(ext1[Result]-ext2[Result]) > TOL) or (Step > TOL) do
 begin
 Step:=Step/2;
 DoStep(y0, y2);
 ext3[Result]:=(ext1[Result]+ext2[Result])/2;
 exitcond(y2, ext3,ExtraData);
 if ext1[Result]*ext3[Result] <= 0 then
  begin
   ext2:=Copy(ext3);
   y1:=Copy(y2);
  end
 else
 begin
  ext1:=Copy(ext3);
   y0:=Copy(y2);
 end;
 end;
 for i:=0 to High(y0) do
  y0[i]:=(y0[i]+y1[i])/2;
  if Assigned(pf) then
  pf(y0, ExtraData);
end;

end.

