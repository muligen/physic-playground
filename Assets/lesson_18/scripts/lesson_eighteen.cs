using System;
using System.Collections.Generic;
using UnityEngine;

enum eField
{
    U_FIELD = 0,
    V_FIELD = 1
}
enum eCellType
{
    FLUID_CELL = 0,
    AIR_CELL = 1,
    SOLID_CELL = 2
}

class FlipFluid
{
    public float density;
    public int fNumX;
    public int fNumY;
    public float h;
    public float fInvSpacing;
    public int fNumCells;
    public float[] u;
    public float[] v;
    public float[] du;
    public float[] dv;
    public float[] prevU;
    public float[] prevV;
    public float[] p;
    public float[] s;
    public eCellType[] cellType;
    public float[] cellColor;
    public int maxParticles;
    public float[] particlePos;
    public float[] particleColor;
    public float[] particleVel;
    public float[] particleDensity;
    public float particleRestDensity;
    public float particleRadius;
    public float pInvSpacing;
    public int pNumX;
    public int pNumY;
    public int pNumCells;
    public int[] numCellParticles;
    public int[] firstCellParticle;
    public int[] cellParticleIds;
    public int numParticles;

    public FlipFluid(float density, float width, float height, float spacing, float particleRadius, int maxParticles)
    {

        // fluid  定义在格子上

        this.density = density;
        this.fNumX = Mathf.FloorToInt(width / spacing) + 1;
        this.fNumY = Mathf.FloorToInt(height / spacing) + 1;
        this.h = Mathf.Max(width / this.fNumX, height / this.fNumY);
        this.fInvSpacing = 1.0f / this.h;
        this.fNumCells = this.fNumX * this.fNumY;

        this.u = new float[this.fNumCells];
        this.v = new float[this.fNumCells];
        this.du = new float[this.fNumCells];
        this.dv = new float[this.fNumCells];
        this.prevU = new float[this.fNumCells];
        this.prevV = new float[this.fNumCells];
        this.p = new float[this.fNumCells];
        this.s = new float[this.fNumCells];
        this.cellType = new eCellType[this.fNumCells];
        this.cellColor = new float[3 * this.fNumCells];

        // particles
        this.maxParticles = maxParticles;

        this.particlePos = new float[2 * this.maxParticles];
        this.particleColor = new float[3 * this.maxParticles];
        for (var i = 0; i < this.maxParticles; i++)
            this.particleColor[3 * i + 2] = 1.0f;

        this.particleVel = new float[2 * this.maxParticles];
        this.particleDensity = new float[this.fNumCells];
        this.particleRestDensity = 0.0f; //密度

        this.particleRadius = particleRadius;
        this.pInvSpacing = 1.0f / (2.2f * particleRadius);
        this.pNumX = Mathf.FloorToInt(width * this.pInvSpacing) + 1;
        this.pNumY = Mathf.FloorToInt(height * this.pInvSpacing) + 1;
        this.pNumCells = this.pNumX * this.pNumY;

        this.numCellParticles = new int[this.pNumCells];
        this.firstCellParticle = new int[this.pNumCells + 1];
        this.cellParticleIds = new int[maxParticles];

        this.numParticles = 0;
    }

    void integrateParticles(float dt, float gravity)
    {
        for (var i = 0; i < this.numParticles; i++)
        {
            this.particleVel[2 * i + 1] += dt * gravity;
            this.particlePos[2 * i] += this.particleVel[2 * i] * dt;
            this.particlePos[2 * i + 1] += this.particleVel[2 * i + 1] * dt;
        }
    }

    void pushParticlesApart(int numIters)
    {
        var colorDiffusionCoeff = 0.001f;   //颜色扩散系数

        // count particles per cell

        Array.Clear(this.numCellParticles, 0, this.numCellParticles.Length);

        for (var i = 0; i < this.numParticles; i++)
        {
            var x = this.particlePos[2 * i];
            var y = this.particlePos[2 * i + 1];

            var xi = Mathf.Clamp(Mathf.FloorToInt(x * this.pInvSpacing), 0, this.pNumX - 1);
            var yi = Mathf.Clamp(Mathf.FloorToInt(y * this.pInvSpacing), 0, this.pNumY - 1);
            var cellNr = xi * this.pNumY + yi;
            this.numCellParticles[cellNr]++;
        }

        // partial sums
        var first = 0;

        for (var i = 0; i < this.pNumCells; i++)
        {
            first += this.numCellParticles[i];
            this.firstCellParticle[i] = first;
        }
        this.firstCellParticle[this.pNumCells] = first;     // guard

        // fill particles into cells

        for (var i = 0; i < this.numParticles; i++)
        {
            var x = this.particlePos[2 * i];
            var y = this.particlePos[2 * i + 1];

            var xi = Mathf.Clamp(Mathf.FloorToInt(x * this.pInvSpacing), 0, this.pNumX - 1);
            var yi = Mathf.Clamp(Mathf.FloorToInt(y * this.pInvSpacing), 0, this.pNumY - 1);
            var cellNr = xi * this.pNumY + yi;
            this.firstCellParticle[cellNr]--;
            this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
        }

        // push particles apart

        var minDist = 2.0f * this.particleRadius;
        var minDist2 = minDist * minDist;

        for (var iter = 0; iter < numIters; iter++)
        {

            for (var i = 0; i < this.numParticles; i++) 
            {
                var px = this.particlePos[2 * i];
                var py = this.particlePos[2 * i + 1];

                var pxi = Mathf.FloorToInt(px * this.pInvSpacing);
                var pyi = Mathf.FloorToInt(py * this.pInvSpacing);
                var x0 = Mathf.Max(pxi - 1, 0);
                var y0 = Mathf.Max(pyi - 1, 0);
                var x1 = Mathf.Min(pxi + 1, this.pNumX - 1);
                var y1 = Mathf.Min(pyi + 1, this.pNumY - 1);

                for (var xi = x0; xi <= x1; xi++)
                {
                    for (var yi = y0; yi <= y1; yi++)
                    {
                        var cellNr = xi * this.pNumY + yi;
                        var begin = this.firstCellParticle[cellNr];
                        var end = this.firstCellParticle[cellNr + 1];
                        for (var j = begin; j < end; j++)
                        {
                            var id = this.cellParticleIds[j];
                            if (id == i)
                                continue;
                            var qx = this.particlePos[2 * id];
                            var qy = this.particlePos[2 * id + 1];

                            var dx = qx - px;
                            var dy = qy - py;
                            var d2 = dx * dx + dy * dy;
                            if (d2 > minDist2 || d2 == 0.0)
                                continue;
                            var d = Mathf.Sqrt(d2);
                            var s = 0.5f * (minDist - d) / d;
                            dx *= s;
                            dy *= s;
                            this.particlePos[2 * i] -= dx;
                            this.particlePos[2 * i + 1] -= dy;
                            this.particlePos[2 * id] += dx;
                            this.particlePos[2 * id + 1] += dy;

                            // diffuse colors   扩散颜色

                            for (var k = 0; k < 3; k++)
                            {
                                var color0 = this.particleColor[3 * i + k];
                                var color1 = this.particleColor[3 * id + k];
                                var color = (color0 + color1) * 0.5f;
                                this.particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
                                this.particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
                            }
                        }
                    }
                }
            }
        }
    }

    void handleParticleCollisions(float obstacleX, float obstacleY, float obstacleRadius)
    {
        var h = 1.0f / this.fInvSpacing;
        var r = this.particleRadius;
        var or = obstacleRadius;
        var or2 = or * or;
        var minDist = obstacleRadius + r;
        var minDist2 = minDist * minDist;

        var minX = h + r;
        var maxX = (this.fNumX - 1) * h - r;
        var minY = h + r;
        var maxY = (this.fNumY - 1) * h - r;


        for (var i = 0; i < this.numParticles; i++)
        {
            var x = this.particlePos[2 * i];
            var y = this.particlePos[2 * i + 1];

            var dx = x - obstacleX;
            var dy = y - obstacleY;
            var d2 = dx * dx + dy * dy;

            // obstacle collision

            if (d2 < minDist2)
            {

                // var d = Mathf.sqrt(d2);
                // var s = (minDist - d) / d;
                // x += dx * s;
                // y += dy * s;

                //this.particleVel[2 * i] = obstacleVelX;
                //this.particleVel[2 * i + 1] = obstacleVelY;  ???
                this.particleVel[2 * i] = 0f;
                this.particleVel[2 * i + 1] = 0f;
            }

            // wall collisions

            if (x < minX)
            {
                x = minX;
                this.particleVel[2 * i] = 0.0f;

            }
            if (x > maxX)
            {
                x = maxX;
                this.particleVel[2 * i] = 0.0f;
            }
            if (y < minY)
            {
                y = minY;
                this.particleVel[2 * i + 1] = 0.0f;
            }
            if (y > maxY)
            {
                y = maxY;
                this.particleVel[2 * i + 1] = 0.0f;
            }
            this.particlePos[2 * i] = x;
            this.particlePos[2 * i + 1] = y;
        }
    }

    void updateParticleDensity()
    {
        var n = this.fNumY;
        var h = this.h;
        var h1 = this.fInvSpacing;
        var h2 = 0.5f * h;

        var d = this.particleDensity; // ??
        Array.Clear(d, 0, d.Length);

        for (var i = 0; i < this.numParticles; i++)
        {
            var x = this.particlePos[2 * i];
            var y = this.particlePos[2 * i + 1];

            x = Mathf.Clamp(x, h, (this.fNumX - 1) * h);
            y = Mathf.Clamp(y, h, (this.fNumY - 1) * h);

            var x0 = Mathf.FloorToInt((x - h2) * h1);
            var tx = ((x - h2) - x0 * h) * h1;
            var x1 = Mathf.Min(x0 + 1, this.fNumX - 2);

            var y0 = Mathf.FloorToInt((y - h2) * h1);
            var ty = ((y - h2) - y0 * h) * h1;
            var y1 = Mathf.Min(y0 + 1, this.fNumY - 2);

            var sx = 1.0f - tx;
            var sy = 1.0f - ty;

            if (x0 < this.fNumX && y0 < this.fNumY) d[x0 * n + y0] += sx * sy;
            if (x1 < this.fNumX && y0 < this.fNumY) d[x1 * n + y0] += tx * sy;
            if (x1 < this.fNumX && y1 < this.fNumY) d[x1 * n + y1] += tx * ty;
            if (x0 < this.fNumX && y1 < this.fNumY) d[x0 * n + y1] += sx * ty;
        }

        if (this.particleRestDensity == 0.0)
        {
            var sum = 0.0f;
            var numFluidCells = 0;

            for (var i = 0; i < this.fNumCells; i++)
            {
                if (this.cellType[i] == eCellType.FLUID_CELL)
                {
                    sum += d[i];
                    numFluidCells++;
                }
            }

            if (numFluidCells > 0)
                this.particleRestDensity = sum / numFluidCells;
        }

        // 			for (var xi = 1; xi < this.fNumX; xi++) {
        // 				for (var yi = 1; yi < this.fNumY; yi++) {
        // 					var cellNr = xi * n + yi;
        // 					if (this.cellType[cellNr] != FLUID_CELL)
        // 						continue;
        // 					var hx = this.h;
        // 					var hy = this.h;

        // 					if (this.cellType[(xi - 1) * n + yi] == SOLID_CELL || this.cellType[(xi + 1) * n + yi] == SOLID_CELL)
        // 						hx -= this.particleRadius;
        // 					if (this.cellType[xi * n + yi - 1] == SOLID_CELL || this.cellType[xi * n + yi + 1] == SOLID_CELL)
        // 						hy -= this.particleRadius;

        // 					var scale = this.h * this.h / (hx * hy)
        // 					d[cellNr] *= scale;
        // 				}
        // 			}
    }

    void transferVelocities(bool toGrid, float flipRatio)
    {
        var n = this.fNumY;
        var h = this.h;
        var h1 = this.fInvSpacing;
        var h2 = 0.5f * h;

        if (toGrid)
        {
            Array.Copy(this.u, this.prevU, this.u.Length);
            Array.Copy(this.v, this.prevV, this.v.Length);

            Array.Clear(this.du, 0, this.du.Length);
            Array.Clear(this.dv, 0, this.du.Length);
            Array.Clear(this.u, 0, this.du.Length);
            Array.Clear(this.v, 0, this.du.Length);

            for (var i = 0; i < this.fNumCells; i++)
                this.cellType[i] = this.s[i] == 0.0f ? eCellType.SOLID_CELL : eCellType.AIR_CELL;

            for (var i = 0; i < this.numParticles; i++)
            {
                var x = this.particlePos[2 * i];
                var y = this.particlePos[2 * i + 1];
                var xi = Mathf.Clamp(Mathf.FloorToInt(x * h1), 0, this.fNumX - 1);
                var yi = Mathf.Clamp(Mathf.FloorToInt(y * h1), 0, this.fNumY - 1);
                var cellNr = xi * n + yi;
                if (this.cellType[cellNr] == eCellType.AIR_CELL)
                    this.cellType[cellNr] = eCellType.FLUID_CELL;
            }
        }

        for (var component = 0; component < 2; component++)
        {

            var dx = component == 0 ? 0.0f : h2;
            var dy = component == 0 ? h2 : 0.0f;

            var f = component == 0 ? this.u : this.v;
            var prevF = component == 0 ? this.prevU : this.prevV;
            var d = component == 0 ? this.du : this.dv;

            for (var i = 0; i < this.numParticles; i++)
            {
                var x = this.particlePos[2 * i];
                var y = this.particlePos[2 * i + 1];

                x = Mathf.Clamp(x, h, (this.fNumX - 1) * h);
                y = Mathf.Clamp(y, h, (this.fNumY - 1) * h);

                var x0 = Mathf.Min(Mathf.FloorToInt((x - dx) * h1), this.fNumX - 2);
                var tx = ((x - dx) - x0 * h) * h1;
                var x1 = Mathf.Min(x0 + 1, this.fNumX - 2);

                var y0 = Mathf.Min(Mathf.FloorToInt((y - dy) * h1), this.fNumY - 2);
                var ty = ((y - dy) - y0 * h) * h1;
                var y1 = Mathf.Min(y0 + 1, this.fNumY - 2);

                var sx = 1.0f - tx;
                var sy = 1.0f - ty;

                var d0 = sx * sy;
                var d1 = tx * sy;
                var d2 = tx * ty;
                var d3 = sx * ty;

                var nr0 = x0 * n + y0;
                var nr1 = x1 * n + y0;
                var nr2 = x1 * n + y1;
                var nr3 = x0 * n + y1;

                if (toGrid)
                {
                    var pv = this.particleVel[2 * i + component];
                    f[nr0] += pv * d0; d[nr0] += d0;
                    f[nr1] += pv * d1; d[nr1] += d1;
                    f[nr2] += pv * d2; d[nr2] += d2;
                    f[nr3] += pv * d3; d[nr3] += d3;
                }
                else
                {
                    var offset = component == 0 ? n : 1;
                    var valid0 = this.cellType[nr0] != eCellType.AIR_CELL || this.cellType[nr0 - offset] != eCellType.AIR_CELL ? 1.0f : 0.0f;
                    var valid1 = this.cellType[nr1] != eCellType.AIR_CELL || this.cellType[nr1 - offset] != eCellType.AIR_CELL ? 1.0f : 0.0f;
                    var valid2 = this.cellType[nr2] != eCellType.AIR_CELL || this.cellType[nr2 - offset] != eCellType.AIR_CELL ? 1.0f : 0.0f;
                    var valid3 = this.cellType[nr3] != eCellType.AIR_CELL || this.cellType[nr3 - offset] != eCellType.AIR_CELL ? 1.0f : 0.0f;

                    var v = this.particleVel[2 * i + component];
                    var dd = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if (dd > 0.0)
                    {

                        var picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / dd;
                        var corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1])
                            + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / dd;
                        var flipV = v + corr;

                        this.particleVel[2 * i + component] = (1.0f - flipRatio) * picV + flipRatio * flipV;
                    }
                }
            }

            if (toGrid)
            {
                for (var i = 0; i < f.Length; i++)
                {
                    if (d[i] > 0.0)
                        f[i] /= d[i];
                }

                // restore solid cells

                for (var i = 0; i < this.fNumX; i++)
                {
                    for (var j = 0; j < this.fNumY; j++)
                    {
                        var solid = this.cellType[i * n + j] == eCellType.SOLID_CELL;
                        if (solid || (i > 0 && this.cellType[(i - 1) * n + j] == eCellType.SOLID_CELL))
                            this.u[i * n + j] = this.prevU[i * n + j];
                        if (solid || (j > 0 && this.cellType[i * n + j - 1] == eCellType.SOLID_CELL))
                            this.v[i * n + j] = this.prevV[i * n + j];
                    }
                }
            }
        }
    }

    void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true)
    {
        Array.Clear(this.p, 0, this.p.Length);
        Array.Copy(this.u, this.prevU, this.prevU.Length);
        Array.Copy(this.v, this.prevV, this.prevV.Length);

        var n = this.fNumY;
        var cp = this.density * this.h / dt;

        for (var i = 0; i < this.fNumCells; i++)
        {
            var u = this.u[i];
            var v = this.v[i];
        }

        for (var iter = 0; iter < numIters; iter++)
        {

            for (var i = 1; i < this.fNumX - 1; i++)
            {
                for (var j = 1; j < this.fNumY - 1; j++)
                {

                    if (this.cellType[i * n + j] != eCellType.FLUID_CELL)
                        continue;

                    var center = i * n + j;
                    var left = (i - 1) * n + j;
                    var right = (i + 1) * n + j;
                    var bottom = i * n + j - 1;
                    var top = i * n + j + 1;

                    var s = this.s[center];
                    var sx0 = this.s[left];
                    var sx1 = this.s[right];
                    var sy0 = this.s[bottom];
                    var sy1 = this.s[top];
                    s = sx0 + sx1 + sy0 + sy1;
                    if (s == 0.0)
                        continue;

                    var div = this.u[right] - this.u[center] +
                        this.v[top] - this.v[center];

                    if (this.particleRestDensity > 0.0 && compensateDrift)
                    {
                        var k = 1.0f;
                        var compression = this.particleDensity[i * n + j] - this.particleRestDensity;
                        if (compression > 0.0)
                            div = div - k * compression;
                    }

                    var p = -div / s;
                    p *= overRelaxation;
                    this.p[center] += cp * p;

                    this.u[center] -= sx0 * p;
                    this.u[right] += sx1 * p;
                    this.v[center] -= sy0 * p;
                    this.v[top] += sy1 * p;
                }
            }
        }
    }

    void updateParticleColors()
    {
        // for (var i = 0; i < this.numParticles; i++) {
        // 	this.particleColor[3 * i] *= 0.99; 
        // 	this.particleColor[3 * i + 1] *= 0.99
        // 	this.particleColor[3 * i + 2] = 
        // 		clamp(this.particleColor[3 * i + 2] + 0.001, 0.0, 1.0)
        // }

        // return;

        var h1 = this.fInvSpacing;

        for (var i = 0; i < this.numParticles; i++)
        {

            var s = 0.01f;

            this.particleColor[3 * i] = Mathf.Clamp(this.particleColor[3 * i] - s, 0.0f, 1.0f);
            this.particleColor[3 * i + 1] = Mathf.Clamp(this.particleColor[3 * i + 1] - s, 0.0f, 1.0f);
            this.particleColor[3 * i + 2] = Mathf.Clamp(this.particleColor[3 * i + 2] + s, 0.0f, 1.0f);

            var x = this.particlePos[2 * i];
            var y = this.particlePos[2 * i + 1];
            var xi = Mathf.Clamp(Mathf.FloorToInt(x * h1), 1, this.fNumX - 1);
            var yi = Mathf.Clamp(Mathf.FloorToInt(y * h1), 1, this.fNumY - 1);
            var cellNr = xi * this.fNumY + yi;

            var d0 = this.particleRestDensity;

            if (d0 > 0.0)
            {
                var relDensity = this.particleDensity[cellNr] / d0;
                if (relDensity < 0.7)
                {
                    s = 0.8f;
                    this.particleColor[3 * i] = s;
                    this.particleColor[3 * i + 1] = s;
                    this.particleColor[3 * i + 2] = 1.0f;
                }
            }
        }
    }

    void setSciColor(int cellNr, float val, float minVal, float maxVal)
    {
        val = Mathf.Min(Mathf.Max(val, minVal), maxVal - 0.0001f);
        var d = maxVal - minVal;
        val = d == 0.0f ? 0.5f : (val - minVal) / d;
        var m = 0.25f;
        var num = Mathf.FloorToInt(val / m);
        var s = (val - num * m) / m;
        float r = 0f, g = 0f, b = 0f;

        switch (num)
        {
            case 0: r = 0.0f; g = s; b = 1.0f; break;
            case 1: r = 0.0f; g = 1.0f; b = 1.0f - s; break;
            case 2: r = s; g = 1.0f; b = 0.0f; break;
            case 3: r = 1.0f; g = 1.0f - s; b = 0.0f; break;
        }

        this.cellColor[3 * cellNr] = r;
        this.cellColor[3 * cellNr + 1] = g;
        this.cellColor[3 * cellNr + 2] = b;
    }

    void updateCellColors()
    {
        Array.Clear(this.cellColor, 0, this.cellColor.Length);

        for (var i = 0; i < this.fNumCells; i++)
        {

            if (this.cellType[i] == eCellType.SOLID_CELL)
            {
                this.cellColor[3 * i] = 0.5f;
                this.cellColor[3 * i + 1] = 0.5f;
                this.cellColor[3 * i + 2] = 0.5f;
            }
            else if (this.cellType[i] == eCellType.FLUID_CELL)
            {
                var d = this.particleDensity[i];
                if (this.particleRestDensity > 0.0)
                    d /= this.particleRestDensity;
                this.setSciColor(i, d, 0.0f, 2.0f);
            }
        }
    }

    public void simulate(float dt, float gravity, float flipRatio, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift, bool separateParticles, float obstacleX, float abstacleY, float obstacleRadius)
    {
        var numSubSteps = 1;
        var sdt = dt / numSubSteps;

        for (var step = 0; step < numSubSteps; step++)
        {
            this.integrateParticles(sdt, gravity);
            if (separateParticles)
                this.pushParticlesApart(numParticleIters);
            this.handleParticleCollisions(obstacleX, abstacleY, obstacleRadius);

            this.transferVelocities(true, flipRatio);
            this.updateParticleDensity();
            this.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
            this.transferVelocities(false, flipRatio);
        }

        this.updateParticleColors();
        this.updateCellColors();

    }
}

class scene
{
    public static float gravity = -9.81f;
    public static float flipRatio = 0.9f;
    public static int numPressureIters = 100;
    public static int numParticleIters = 2;
    public static float overRelaxation = 1.9f;
    public static bool compensateDrift = true;
    public static bool separateParticles = true;
    public static float obstacleX = 0.0f;
    public static float obstacleY = 0.0f;
    public static float obstacleRadius= 0.15f;
    public static bool paused = true;
    public static bool showObstacle = true;
    public static float obstacleVelX= 0.0f;
    public static float obstacleVelY= 0.0f;
    public static bool showParticles = true;
    public static bool showGrid= false;
}


public class lesson_eighteen : MonoBehaviour
{
    FlipFluid f;
    public GameObject prefab;
    public GameObject obstacle;
    List<GameObject> paticles;

    void setObstacle(float x, float y, FlipFluid f, bool reset )
    {
        var vx = 0.0f;
        var vy = 0.0f;
        if (!reset)
        {
            vx = (x - scene.obstacleX) / Time.fixedDeltaTime;
            vy = (y - scene.obstacleY) / Time.fixedDeltaTime;
        }
        scene.obstacleX = x;
        scene.obstacleY = y;
        var r = scene.obstacleRadius;
        var n = f.fNumY;
        var cd = Mathf.Sqrt(2) * f.h;

        for (var i = 1; i < f.fNumX - 2; i++)
        {
            for (var j = 1; j < f.fNumY - 2; j++)
            {

                f.s[i * n + j] = 1.0f;

                var dx = (i + 0.5) * f.h - x;
                var dy = (j + 0.5) * f.h - y;

                if (dx * dx + dy * dy < r * r)
                {
                    f.s[i * n + j] = 0.0f;
                    f.u[i * n + j] = vx;
                    f.u[(i + 1) * n + j] = vx;
                    f.v[i * n + j] = vy;
                    f.v[i * n + j + 1] = vy;
                }
            }
        }

        scene.showObstacle = true;
        scene.obstacleVelX = vx;
        scene.obstacleVelY = vy;
    }

    void Start()
    {
        scene.obstacleRadius = 0.15f;
        scene.overRelaxation = 1.9f;

        scene.numPressureIters = 50;
        scene.numParticleIters = 2;

        var res = 100f;

        var tankHeight = 2.5f;
        var tankWidth = 4f;
        var h = tankHeight / res;
        var density = 1000.0f;

        var relWaterHeight = 0.3f;

        var relWaterWidth = 0.3f;

        // dam break

        // compute number of particles

        var r = 0.3f * h;    // particle radius w.r.t. cell size
        var dx = 2.0f * r;
        var dy = Mathf.Sqrt(3.0f) / 2.0f * dx;

        var numX = Mathf.FloorToInt((relWaterWidth * tankWidth - 2.0f * h - 2.0f * r) / dx);
        var numY = Mathf.FloorToInt((relWaterHeight * tankHeight - 2.0f * h - 2.0f * r) / dy);
        var maxParticles = numX * numY;

        // create fluid

        f = new FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

        // create particles
 
        f.numParticles = numX * numY;
        var p = 0;
        for (var i = 0; i < numX; i++)
        {
            for (var j = 0; j < numY; j++)
            {
                f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0f : r);
                f.particlePos[p++] = h + r + dy * j;
            }
        }

        // setup grid cells for tank

        var n = f.fNumY;

        for (var i = 0; i < f.fNumX; i++)
        {
            for (var j = 0; j < f.fNumY; j++)
            {
                var s = 1.0f;    // fluid
                if (i == 0 || i == f.fNumX - 1 || j == 0)
                    s = 0.0f;    // solid
                f.s[i * n + j] = s;

            }
        }

        setObstacle(2.0f, 1.0f, f, true);
        paticles = new List<GameObject>();
        for (int i = 0; i < f.numParticles; i++)
        {
            paticles.Add(Instantiate(prefab, new Vector3(f.particlePos[i * 2], f.particlePos[i * 2 + 1], 0), Quaternion.identity));
        }
        Debug.Log(paticles.Count);
    }

    private void FixedUpdate()
    {
        f.simulate(
            Time.fixedDeltaTime, scene.gravity, scene.flipRatio, scene.numPressureIters, scene.numParticleIters,
            scene.overRelaxation, scene.compensateDrift, scene.separateParticles,
            scene.obstacleX, scene.obstacleY, scene.obstacleRadius);

        //draw
        for (int i = 0; i < paticles.Count; i++)
        {
            paticles[i].GetComponent<Transform>().position = new Vector3(f.particlePos[i * 2], f.particlePos[i * 2 + 1], 0);
            paticles[i].GetComponent<MeshRenderer>().material.color = new Color(f.particleColor[i * 3], f.particleColor[i * 3 + 1], f.particleColor[i * 3 + 2]);
        }

        if (Input.GetMouseButton(0))
        {
            setObstacle(Input.mousePosition.x / 400f, Input.mousePosition.y / 400f, f, false);
        }
        obstacle.GetComponent<Transform>().position = new Vector3(scene.obstacleX, scene.obstacleY, 0);
    }

}
