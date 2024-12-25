class OptimumResult
{
    public Vector xopt { get; set; }
    public Vector? ogr { get; set; }
    public Vector fopt { get; set; }
    public int CountIterations = 0;
    public OptimumResult(Vector x, Vector? ogr, Vector f, int CountIter = 0)
    {
        xopt = x.Copy();
        this.ogr = null;
        if (ogr != null) this.ogr = ogr.Copy();
        fopt = f?.Copy();
        CountIterations = CountIter;
    }

    public override string ToString()
    {
        if (ogr == null) return string.Format(" Result xopt={0}  fopt={1} count {2}", xopt, fopt, CountIterations);
        return string.Format(" Result xopt={0}  ogr={1}  fopt={2} count {3}", xopt, ogr, fopt, CountIterations);
    }
}