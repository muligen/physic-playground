using System;
using UnityEngine;
using UnityEngine.UI;

class Fluid
{
    public float density;
    public int numX;
    public int numY;
    public int numCells;
    public float h;
    public float[] u;
    public float[] v;
    public float[] newU;
    public float[] newV;
    public float[] p;
    public float[] s;
    public float[] m;
    public float[] newM;

    public Fluid(float density, int numX, int numY, float h)
    {
        this.density = density;
        this.numX = numX + 2;
        this.numY = numY + 2;
        this.numCells = this.numX * this.numY;
        this.h = h;
        this.u = new float[this.numCells];
        this.v = new float[this.numCells];
        this.newU = new float[this.numCells];
        this.newV = new float[this.numCells];
        this.p = new float[this.numCells];
        this.s = new float[this.numCells];
        this.m = new float[this.numCells];
        this.newM = new float[this.numCells];
        for (int i = 0; i < this.m.Length; i++)
        {
            this.m[i] = 1.0f;
        }

        var num = numX * numY;
    }

    public void integrate(float dt, float gravity)
    {
        var n = this.numY;
        for (var i = 1; i < this.numX; i++)
        {
            for (var j = 1; j < this.numY - 1; j++)
            {
                if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0)
                    this.v[i * n + j] += gravity * dt;
            }
        }
    }

    public void solveIncompressibility(int numIters, float dt)
    {

        var n = this.numY;
        var cp = this.density * this.h / dt;

        for (var iter = 0; iter < numIters; iter++)
        {

            for (var i = 1; i < this.numX - 1; i++)
            {
                for (var j = 1; j < this.numY - 1; j++)
                {

                    if (this.s[i * n + j] == 0.0)
                        continue;

                    var s = this.s[i * n + j];
                    var sx0 = this.s[(i - 1) * n + j];
                    var sx1 = this.s[(i + 1) * n + j];
                    var sy0 = this.s[i * n + j - 1];
                    var sy1 = this.s[i * n + j + 1];
                    s = sx0 + sx1 + sy0 + sy1;
                    if (s == 0.0)
                        continue;

                    var div = this.u[(i + 1) * n + j] - this.u[i * n + j] +
                        this.v[i * n + j + 1] - this.v[i * n + j];

                    var p = -div / s;
                    p *= 1.9f; //  1到2，可以加速收敛的参数。
                    this.p[i * n + j] += cp * p;

                    this.u[i * n + j] -= sx0 * p;
                    this.u[(i + 1) * n + j] += sx1 * p;
                    this.v[i * n + j] -= sy0 * p;
                    this.v[i * n + j + 1] += sy1 * p;
                }
            }
        }
    }

    public void extrapolate() // 假定边界速度始终持平
    {
        var n = this.numY;
        for (var i = 0; i < this.numX; i++)
        {
            this.u[i * n + 0] = this.u[i * n + 1];
            this.u[i * n + this.numY - 1] = this.u[i * n + this.numY - 2];
        }
        for (var j = 0; j < this.numY; j++)
        {
            this.v[0 * n + j] = this.v[1 * n + j];
            this.v[(this.numX - 1) * n + j] = this.v[(this.numX - 2) * n + j];
        }
    }

    public float sampleField(float x, float y, int field)
    {
        var n = this.numY;
        var h = this.h;
        var h1 = 1.0f / h;
        var h2 = 0.5f * h;

        x = Mathf.Max(Mathf.Min(x, this.numX * h), h);
        y = Mathf.Max(Mathf.Min(y, this.numY * h), h);

        var dx = 0.0f;
        var dy = 0.0f;

        float[] f = {};

        switch (field)
        {
            case 0: f = this.u; dy = h2; break;
            case 1: f = this.v; dx = h2; break;
            case 2:
                f = this.m; dx = h2; dy = h2; break;
        }

        var x0 = Mathf.Min(Mathf.FloorToInt((x - dx) * h1), this.numX - 1);
        var tx = ((x - dx) - x0 * h) * h1;
        var x1 = Mathf.Min(x0 + 1, this.numX - 1);

        var y0 = Mathf.Min(Mathf.FloorToInt((y - dy) * h1), this.numY - 1);
        var ty = ((y - dy) - y0 * h) * h1;
        var y1 = Mathf.Min(y0 + 1, this.numY - 1);

        var sx = 1.0f - tx;
        var sy = 1.0f - ty;

        var val = sx * sy * f[x0 * n + y0] +
                    tx * sy * f[x1 * n + y0] +
                    tx * ty * f[x1 * n + y1] +
                    sx * ty * f[x0 * n + y1];

        return val;
    }

    public float avgU(int i, int j)
    {
        var n = this.numY;
        var u = (this.u[i * n + j - 1] + this.u[i * n + j] +
            this.u[(i + 1) * n + j - 1] + this.u[(i + 1) * n + j]) * 0.25f;
        return u;
    }

    public float avgV(int i, int j)
    {
        var n = this.numY;
        var v = (this.v[(i - 1) * n + j] + this.v[i * n + j] +
            this.v[(i - 1) * n + j + 1] + this.v[i * n + j + 1]) * 0.25f;
        return v;
    }

    public void advectVel(float dt)
    {
        for (int i = 0; i < this.numCells; i++)
        {
            this.newU[i] = this.u[i];
            this.newV[i] = this.v[i];
        }

        var n = this.numY;
        var h = this.h;
        var h2 = 0.5f * h;

        for (var i = 1; i < this.numX; i++)
        {
            for (var j = 1; j < this.numY; j++)
            {
                // u component
                if (this.s[i * n + j] != 0.0 && this.s[(i - 1) * n + j] != 0.0 && j < this.numY - 1)
                {
                    var x = i * h;
                    var y = j * h + h2;
                    var u = this.u[i * n + j];
                    var v = this.avgV(i, j);

                    x = x - dt * u;
                    y = y - dt * v;
                    u = this.sampleField(x, y, 0);
                    this.newU[i * n + j] = u;
                }
                // v component
                if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0 && i < this.numX - 1)
                {
                    var x = i * h + h2;
                    var y = j * h;
                    var u = this.avgU(i, j);

                    var v = this.v[i * n + j];
                    x = x - dt * u;
                    y = y - dt * v;
                    v = this.sampleField(x, y, 1);
                    this.newV[i * n + j] = v;
                }
            }
        }
        for (int i = 0; i < this.numCells; i++)
        {
            this.u[i] = this.newU[i];
            this.v[i] = this.newV[i];
        }
    }

    public void advectSmoke(float dt)
    {
        for (int i = 0; i < this.numCells; i++)
        {
            this.newM[i] = this.m[i];
        }

        var n = this.numY;
        var h = this.h;
        var h2 = 0.5f * h;

        for (var i = 1; i < this.numX - 1; i++)
        {
            for (var j = 1; j < this.numY - 1; j++)
            {

                if (this.s[i * n + j] != 0.0)
                {
                    var u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5f;
                    var v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5f;
                    var x = i * h + h2 - dt * u;
                    var y = j * h + h2 - dt * v;

                    this.newM[i * n + j] = this.sampleField(x, y, 2);
                }
            }
        }
        for (int i = 0; i < this.numCells; i++)
        {
            this.m[i] = this.newM[i];
        }
    }


    public void simulate(float dt, float gravity, int numIters)
    {

        this.integrate(dt, gravity);

        Array.Clear(this.p, 0, this.p.Length);
        this.solveIncompressibility(numIters, dt);

        this.extrapolate();
        this.advectVel(dt);
        this.advectSmoke(dt);
    }
}

public class lesson_seventeen : MonoBehaviour
{
    public Image img_Fluid2D;
    public Image obstacle;
    int numX;
    int numY;
    Fluid f;
    void Start()
    {
        numX = img_Fluid2D.GetComponent<drawFluid2d>().numX;
        numY = img_Fluid2D.GetComponent<drawFluid2d>().numY;
        f = new Fluid(1.0f, numX-2, numY-2, 1.0f/numY);
        img_Fluid2D.GetComponent<drawFluid2d>().smoke_c = f.m;  // 传一个引用过去
        SetUpScene();
    }
    void setObstacle(Fluid f)
    {

        var x = (obstacle.transform.position.x - img_Fluid2D.transform.position.x) / (numY * 10);
        var y = (obstacle.transform.position.y - img_Fluid2D.transform.position.y + 10f) / (numY * 10);
        
        var vx = 0.0f;
        var vy = 0.0f;

        var r = 10f / numY;
        var n = f.numY;

        for (var i = 1; i < f.numX - 2; i++)
        {
            for (var j = 1; j < f.numY - 2; j++)
            {

                f.s[i * n + j] = 1.0f;

                var dx = (i + 0.5f) * f.h - x;
                var dy = (j + 0.5f) * f.h - y;

                if (dx * dx + dy * dy < r * r)
                {
                    f.s[i * n + j] = 0.0f;
                    f.m[i * n + j] = 1.0f;
                    f.u[i * n + j] = vx;
                    f.u[(i + 1) * n + j] = vx;
                    f.v[i * n + j] = vy;
                    f.v[i * n + j + 1] = vy;
                }
            }
        }
    }
    void SetUpScene()
    {
        int n = numY;
        var inVel = 2.0f;
        for (var i = 0; i < f.numX; i++)
        {
            for (var j = 0; j < f.numY; j++)
            {
                var s = 1.0f;    // fluid
                if (i == 0 || j == 0 || j == f.numY - 1)
                    s = 0.0f;    // solid
                f.s[i * n + j] = s;
                if (i == 1)
                {
                    f.u[i * n + j] = inVel;
                }
            } 
        }

        var pipeH = 0.1f * f.numY;
        var minJ = Mathf.FloorToInt(0.5f * f.numY - 0.5f * pipeH);
        var maxJ = Mathf.FloorToInt(0.5f * f.numY + 0.5f * pipeH);

        for (var j = minJ; j < maxJ; j++)
            f.m[j] = 0.0f;

        setObstacle(f);

    }

    private void FixedUpdate()
    {
        f.simulate(Time.fixedDeltaTime, 0f, 100);
        img_Fluid2D.GetComponent<drawFluid2d>().SetVerticesDirty();
        if (Input.GetMouseButton(0))
        {
            obstacle.transform.position = Input.mousePosition;
            setObstacle(f);
        }
    }
}
