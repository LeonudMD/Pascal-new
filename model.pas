{$mode ObjFPC}{$H+}
unit Model;
interface
uses uSolver;
const muTerra = 398600.44158;  // гравитационный параметр Земли, км^3/с^2
      Earth_C20 = -1098.08E-6; // коэффициент при второй зональной гармонике
                               // в разложении потенциала по сферическим функциям
      Earth_eqR = 6378.16;     // экваториальный радиус Земли, км

      // звёздные сутки, сек
      Earth_Star_Period = 23*3600 + 56*60 + 04.09054;

//параметры задачи
type TPerturbation = (ptGrav20,      // возмущения от второй зональной гармоники, км/с
                      ptPropulsion   // возмущения от двигателя, км/c
                      );

     TTaskParam = record
      Perturbation:  set of TPerturbation; // множество учитываемых возмущений
      begin_accel: real;  // начальное реактивное ускорение, км/c^2
      rate_outflow: real; // скорость истечения рабочего тела, км/с

      dinclination: real; // разность наклонений.(ik - i0)
      dradius: real;      // отношение радиусов орбит. (rk/r0)
      Vxk: real;          // характеристическая скорость перелета, км/с
     end;

     //
     TMTAParam = record
      AlphaEY,   // удельная масса энергоустановки, кг/кВт
      NuT,       // тяговый кпд ДУ, %
      NuEY,      // кпд преобразователя энергии, %
      GammaDY,   // удельная масса ДУ, кг/Н
      GammaSPX,  // удельная масса системы хранения и подачи раб. тела
      GammaK   :real; // относительная масса конструкции
     end;

{ возвращает мощность энергоустановки, кг
  }
function GetMTAPower(const p: TMTAParam; rate, thrust: real): real;

{ возвращает стартовую массу МТА, кг
 }
function GetStartMass(const p: TMTAParam; rate, Tm,Vxk, mpn: real): real;

{ возвращает оптимальную скорость истечения рабочего тела, км/с
}
function GetRateopt(const p: TMTAParam; Tm:real):real;

{возвращает характеристическую скорость перелёта, км/с
r0 — радиус начальной орбиты, км
i0 — наклонение начальной орбиты, рад
rk — радиус конечной орбиты, км
ik — наклонение конечной орбиты, км }
function GetCharacteristicSpeed(r0, i0, rk, ik: real): real;

{инициализация параметров задачи
r0 — радиус начальной орбиты, км
i0 — наклонение начальной орбиты, рад
rk — радиус конечной орбиты, км
ik — наклонение конечной орбиты, км
Tm — время перелёта, сек
с  — скорость истечения рабочего тела, км/с

ВОЗВРАЩАЕМОЕ ЗНАЧЕНИЕ
Task — записть с расчитанными параметрами задачи. }
procedure InitTaskParam(out Task: TTaskParam; r0, i0, rk, ik, Tm, c: real);

{возвращает период обращения, сек
 a - большая полуось, км}
function GetOrbitalPeriod(a: real): real;

{возвращает большую полуось орбиты по периоду обращения, км
period — орбитальный период обращения, сек}
function GetOrbitBigSemiAxe(period: real): real;

{возвращает значение амплитуды угла управления
 Task — параметры задачи. Должны быть предварительно рассчитаны с помощью
        InitTaskParam, например.
 Vx   — характеристическая скорость перелёта, км/c }
function GetPsim(const Task: TTaskParam; Vx: real): real;

{возвращает текущую характеристическую скорость
 Task — параметры задачи. Должны быть предварительно рассчитаны с помощью
        InitTaskParam, например.
 t — текущее время, с
}
function GetCurrentChSpeed(const Task: TTaskParam; t: real): real;

{условие прерывания интегрирования}
procedure ExitCondition(const inv: TVector; var outv: TVector;
                          ExtraData: pointer = nil);

{правая часть ОДУ модели}
procedure ModelRp(const inv: TVector; var outv: TVector;
                          ExtraData: pointer = nil);

implementation
uses Math;

function GetMTAPower(const p: TMTAParam; rate, thrust: real): real;
begin
  Result:=thrust*rate/(2*p.NuT*p.NuEY);
end;



function GetRateopt(const p: TMTAParam; Tm: real): real;
begin
  Result:=sqrt(Tm*p.NuT*p.NuEY*(1+p.GammaSPX)/(p.AlphaEY*0.001+p.GammaK));
end;

function GetCharacteristicSpeed(r0, i0, rk, ik: real): real;
var dk, dr: real;
begin
 dk:=0.5*pi*(ik - i0);
 dr:=sqrt(rk/r0);
 Result:=sqrt(muTerra/r0)*sqrt(1 -2*cos(dk)/dr + r0/rk);
end;

function GetStartMass(const p: TMTAParam; rate, Tm, Vxk, mpn: real): real;
var z1, Mu, A: real;
begin
  z1:=1-exp(-Vxk/rate);
  A:=1+p.GammaSPX-(rate)/Tm*(p.GammaDY+p.GammaK+(rate)/(2*p.NuT*p.NuEY)*p.AlphaEY*0.001);
  Mu:=1-A*z1;
 Result:=mpn/Mu;
end;

procedure InitTaskParam(out Task: TTaskParam; r0, i0, rk, ik, Tm, c: real);
begin
 Task.Perturbation:=[];
 Task.Vxk:=GetCharacteristicSpeed(r0, i0, rk, ik);
 Task.begin_accel:=c/Tm*(1 - exp(-Task.Vxk/c));
 Task.dinclination:=ik - i0;
 Task.dradius:=rk/r0;
 Task.rate_outflow:=c;
end;

function GetOrbitalPeriod(a: real): real;
begin
 Result:=2*pi*sqrt(a*a*a/muTerra);
end;

function GetOrbitBigSemiAxe(period: real): real;
begin
 Result:=Power(muTerra*sqr(0.5*period/pi), 1/3);
end;

function GetPsim(const Task: TTaskParam; Vx: real): real;
{
 psim = atan(a/(1 - b - Vx*c^2))
}
const eps: real = 1E-15; // точность детектирования нуля
var a, b, c: real; // параметры программы управления
    arg: real;     // значение знаменателя
begin
 a:=sin(0.5*Task.dinclination*pi)/sqrt(Task.dradius);
 b:=cos(0.5*Task.dinclination*pi)/sqrt(Task.dradius);
 c:=sqrt(1 - 2*b + 1/Task.dradius);

 arg:=1 - b - Vx/Task.Vxk*sqr(c);
 if abs(arg) < eps then   // если в знаменателе почти ноль, psim = pi/2
  Result:=pi/2
 else if arg > 0 then
  Result:=-arctan(a/arg)
 else
  Result:=pi - arctan(a/arg);
end;

function GetCurrentChSpeed(const Task: TTaskParam; t: real): real;
begin
 with Task do
  Result:=-rate_outflow*ln(1 - begin_accel/rate_outflow*t);
end;

procedure ExitCondition(const inv: TVector; var outv: TVector; ExtraData: pointer);
var Vx: real; // текущая характеристическая скорость, км/c
begin
 Assert(ExtraData <> nil, 'Не заданы параметры задачи!');

 {прерывание интегрирования по характеристической скорости}
 Vx:=GetCurrentChSpeed(TTaskParam(ExtraData^), inv[0]);
 outv[0]:= Vx - TTaskParam(ExtraData^).Vxk;
end;

procedure ModelRp(const inv: TVector; var outv: TVector; ExtraData: pointer);
      // коэффициент для расчёта возмущений от второй зональной гармоники км^5/с^2
const Earth_eps = -0.5*3*Earth_C20*sqr(Earth_eqR)*muTerra;
var W, S, T: real; // проекции возмущаещего ускорения, км/c
    r:  real;      // радиус вектор, км

    cosu, sinu,
    cost, sint,
    cosi, sini,
    cos_psi, sin_psi: real;

    tmp, tmp2: real;
    Vx: real;      // текущая характеристическая скорость
    ap: real;      // ускорение от двигательной установки, км/с
    psim: real;    // амплитуда угла управления, рад
begin
 Assert(ExtraData <> nil, 'Не заданы параметры задачи!');

 {0   1   2 3   4   5   6
  t Omega i p omega e tetta}
 W:=0; S:=0;  T:=0;

 // вычисление sin и cos углов
 SinCos(inv[4] + inv[6], sinu, cosu);  // u = omega + tetta
 SinCos(inv[6], sint, cost);           // tetta
 SinCos(inv[2], sini, cosi);           // inclination

 r:=inv[3]/(1 + inv[5]*cost);          // r — радиус-вектор
 tmp:=sqrt(inv[3]/muTerra);
 tmp2:=r/sqrt(muTerra*inv[3]);

 // вычисление возмущающих ускорений
 // возмущающее ускорение от второй зональной гармоники
 if ptGrav20 in TTaskParam(ExtraData^).Perturbation then begin
  S:=S + Earth_eps/sqr(sqr(r))*(3*sqr(sini)*sqr(sinu) - 1);
  T:=T - 2*Earth_eps/sqr(sqr(r))*sqr(sini)*sinu*cosu;
  W:=W - 2*Earth_eps/sqr(sqr(r))*sini*cosi*sinu;
 end;

 // возмущение от двигательной установки
 if ptPropulsion in TTaskParam(ExtraData^).Perturbation then begin
  // ускорение от двигательной установки, км/с
  with TTaskParam(ExtraData^) do
   ap:=begin_accel * rate_outflow / (rate_outflow - begin_accel*inv[0]);

  // текущая характеристическая скорость
  Vx:=GetCurrentChSpeed(TTaskParam(ExtraData^), inv[0]);

  // амплитуда угла управления
  psim:=GetPsim(TTaskParam(ExtraData^), Vx);
  SinCos(psim, sin_psi, cos_psi);

  T:=T + ap*cos_psi;
  W:=W + ap*sin_psi*sign(-cosu);
 end;

 // вычисление правых частей ОДУ
 outv[1]:=tmp2*W*sinu/sini;
 outv[2]:=tmp2*W*cosu;
 outv[3]:=tmp*2*T*r;
 outv[4]:=tmp*(-S*cost/inv[5] + T*(1 + r/inv[3])*sint*T/inv[5] -
                W*r/inv[3]*cosi/sini*sinu);
 outv[5]:=tmp*(S*sint + T*(inv[5]*r/inv[3] + (1 + r/inv[3])*cost));
 outv[6]:=tmp*(muTerra/sqr(r) + S*cost/inv[5] - T*(1 + r/inv[3])*sint/inv[5]);
end;

end.


