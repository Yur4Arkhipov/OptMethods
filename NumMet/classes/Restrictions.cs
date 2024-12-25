
enum TypeRestrictions
{
    Argument,
    Functions,
    ArgumentsAndFunctions
}
public enum RestrictionTypes
{
    DoubleOgr,
    UpperOgr,
    LowerOgr
}
// public delegate bool DelTestOgr(Vector x);


// Простое скалярное ограничение
class Restriction
{
    public double Upper { get; set; }
    public double Lower { get; set; }
    public RestrictionTypes Type
    {
        get; set;
    }

    public Restriction(double upper, double lower, RestrictionTypes type)
    {
        Upper = upper;
        Lower = lower;
        Type = type;
    }
    public Restriction()
    {
        Upper = double.MaxValue;
        Lower = double.MinValue;
        Type = RestrictionTypes.DoubleOgr;
    }
    public bool IsRegion(double x)
    {
        if (Type == RestrictionTypes.UpperOgr)
        {
            if (x <= Upper) return true;
        }
        if (Type == RestrictionTypes.LowerOgr)
        {
            if (x >= Lower) return true;
        }
        if (Type == RestrictionTypes.DoubleOgr)
        {
            if (x >= Lower && x <= Upper) return true;
        }
        return false;
    }
    public double Proection(double x)
    {
        double xp = x;
        if (Type == RestrictionTypes.UpperOgr)
        {
            if (x > Upper) xp = Upper;
        }
        if (Type == RestrictionTypes.LowerOgr)
        {
            if (x < Lower) xp = Lower;
        }
        if (Type == RestrictionTypes.DoubleOgr)
        {
            if (x < Lower) xp = Lower;
            if (x > Upper) xp = Upper;
        }
        return xp;
    }
}


// Функциональное ограничение
class RestrictionFunc
{
    public double Upper { get; set; }
    public double Lower { get; set; }
    Func<Vector, double> function;
    public RestrictionTypes Type
    {
        get; set;
    }

    public RestrictionFunc(double upper, double lower, RestrictionTypes type, Func<Vector, double> func)
    {
        Upper = upper;
        Lower = lower;
        Type = type;
        function = func;
    }
    public RestrictionFunc(Func<Vector, double> func)
    {
        Upper = double.MaxValue;
        Lower = double.MinValue;
        Type = RestrictionTypes.DoubleOgr;
        function = func;
    }
    public bool IsRegion(Vector v)
    {
        double x = function(v);
        if (Type == RestrictionTypes.UpperOgr)
        {
            if (x <= Upper) return true;
        }
        if (Type == RestrictionTypes.LowerOgr)
        {
            if (x >= Lower) return true;
        }
        if (Type == RestrictionTypes.DoubleOgr)
        {
            if (x >= Lower && x <= Upper) return true;
        }
        return false;
    }

    public double GetValue(Vector v)
    {
        return function(v);
    }
}