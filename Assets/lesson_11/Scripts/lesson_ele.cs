using System;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

class Hash
{
    public float spacing;
    public int tableSize;
    public int[] cellStart;
    public int[] cellEntries;
    public int[] queryIds;
    public int querySize;

    public Hash(float spacing, int maxNumObjects)
    {
        this.spacing = spacing;
        this.tableSize = 2 * maxNumObjects;
        this.cellStart = new int[this.tableSize + 1];
        this.cellEntries = new int[maxNumObjects];
        this.queryIds = new int[maxNumObjects];
        this.querySize = 0;
    }

    public int hashCoords(int xi, int yi, int zi)
    {
        var h = (xi * 92837111) ^ (yi * 689287499) ^ (zi * 283923481);  // fantasy function
        return Mathf.Abs(h) % this.tableSize;
    }

    public int intCoord(float coord)
    {
        return Mathf.FloorToInt(Mathf.Floor(coord / this.spacing));
    }

    public int hashPos(float[] pos, int nr)
    {
        return this.hashCoords(
            this.intCoord(pos[3 * nr]),
            this.intCoord(pos[3 * nr + 1]),
        this.intCoord(pos[3 * nr + 2]));
    }

    public void create(float[] pos)
    {
        var numObjects = Mathf.Min(pos.Length / 3, this.cellEntries.Length);

        // determine cell sizes

        Array.Clear(this.cellStart, 0, this.cellStart.Length);
        //this.cellStart.fill(0);
        //this.cellEntries.fill(0);

        for (var i = 0; i < numObjects; i++)
        {
            var h = this.hashPos(pos, i);
            this.cellStart[h]++;
        }

        // determine cells starts

        var start = 0;
        for (var i = 0; i < this.tableSize; i++)
        {
            start += this.cellStart[i];
            this.cellStart[i] = start;
        }
        this.cellStart[this.tableSize] = start; // guard

        // fill in objects ids

        for (var i = 0; i < numObjects; i++)
        {
            var h = this.hashPos(pos, i);
            this.cellStart[h]--;
            this.cellEntries[this.cellStart[h]] = i;
        }
    }

    public void query(float[] pos, int nr, float maxDist)
    {
        var x0 = this.intCoord(pos[3 * nr] - maxDist);
        var y0 = this.intCoord(pos[3 * nr + 1] - maxDist);
        var z0 = this.intCoord(pos[3 * nr + 2] - maxDist);

        var x1 = this.intCoord(pos[3 * nr] + maxDist);
        var y1 = this.intCoord(pos[3 * nr + 1] + maxDist);
        var z1 = this.intCoord(pos[3 * nr + 2] + maxDist);

        this.querySize = 0;

        for (var xi = x0; xi <= x1; xi++)
        {
            for (var yi = y0; yi <= y1; yi++)
            {
                for (var zi = z0; zi <= z1; zi++)
                {
                    var h = this.hashCoords(xi, yi, zi);
                    var start = this.cellStart[h];
                    var end = this.cellStart[h + 1];

                    for (var i = start; i < end; i++)
                    {
                        this.queryIds[this.querySize] = this.cellEntries[i];
                        this.querySize++;
                    }
                }
            }
        }
    }
}class ThreeDBalls
{
    private float radius;
    private float[] pos;
    private float[] prevPos;
    private float[] vel;
    private int numBalls;
    private Hash hash;
    private float[] normal;
    private List<GameObject> balls;

    #region 向量计算辅助方法
    void vecSetZero(float[] a, int anr)
    {
        anr *= 3;
        a[anr++] = 0.0f;
        a[anr++] = 0.0f;
        a[anr] = 0.0f;
    }

    void vecScale(float[] a, int anr, float scale)
    {
        anr *= 3;
        a[anr++] *= scale;
        a[anr++] *= scale;
        a[anr] *= scale;
    }

    void vecCopy(float[] a, int anr, float[] b, int bnr)
    {
        anr *= 3; bnr *= 3;
        a[anr++] = b[bnr++];
        a[anr++] = b[bnr++];
        a[anr] = b[bnr];
    }

    void vecAdd(float[] a, int anr, float[] b, int bnr, float scale = 1.0f)
    {
        anr *= 3; bnr *= 3;
        a[anr++] += b[bnr++] * scale;
        a[anr++] += b[bnr++] * scale;
        a[anr] += b[bnr] * scale;
    }

    void vecSetDiff(float[] dst, int dnr, float[] a, int anr, float[] b, int bnr, float scale = 1.0f)
    {
        dnr *= 3; anr *= 3; bnr *= 3;
        dst[dnr++] = (a[anr++] - b[bnr++]) * scale;
        dst[dnr++] = (a[anr++] - b[bnr++]) * scale;
        dst[dnr] = (a[anr] - b[bnr]) * scale;
    }

    float vecLengthSquared(float[] a, int anr)
    {
        anr *= 3;
        float a0 = a[anr], a1 = a[anr + 1], a2 = a[anr + 2];
        return a0 * a0 + a1 * a1 + a2 * a2;
    }

    float vecDistSquared(float[] a, int anr, float[] b, int bnr)
    {
        anr *= 3; bnr *= 3;
        float a0 = a[anr] - b[bnr], a1 = a[anr + 1] - b[bnr + 1], a2 = a[anr + 2] - b[bnr + 2];
        return a0 * a0 + a1 * a1 + a2 * a2;
    }

    float vecDot(float[] a, int anr, float[] b, int bnr)
    {
        anr *= 3; bnr *= 3;
        return a[anr] * b[bnr] + a[anr + 1] * b[bnr + 1] + a[anr + 2] * b[bnr + 2];
    }

    void vecSetCross(float[] a, int anr, float[] b, int bnr, float[] c, int cnr)
    {
        anr *= 3; bnr *= 3; cnr *= 3;
        a[anr++] = b[bnr + 1] * c[cnr + 2] - b[bnr + 2] * c[cnr + 1];
        a[anr++] = b[bnr + 2] * c[cnr + 0] - b[bnr + 0] * c[cnr + 2];
        a[anr] = b[bnr + 0] * c[cnr + 1] - b[bnr + 1] * c[cnr + 0];
    }
    #endregion
    public ThreeDBalls(float radius, float[] pos, float[] vel)
    {
        // physics data 

        this.radius = radius;
        this.pos = pos;
        this.prevPos = pos;
        this.vel = vel;
        this.numBalls = pos.Length / 3;
        this.hash = new Hash(2.0f * radius, this.numBalls);
        //this.showCollisions = false;

        this.normal = new float[3];

        var mat = Resources.Load<Material>("bunny/mat_pink");
        this.balls = new List<GameObject>();

        for (int i = 0; i < this.pos.Length/3; i++)
        {
            var surfaceMesh = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            surfaceMesh.GetComponent<MeshRenderer>().material = mat;
            surfaceMesh.transform.position = new Vector3(this.pos[i * 3], this.pos[i * 3 + 1], this.pos[i * 3 + 2]);
            surfaceMesh.transform.localScale = new Vector3(this.radius, this.radius, this.radius);
            this.balls.Add(surfaceMesh);
        }
    }

    private void updateMesh()
    {
        for (int i = 0; i < this.balls.Count; i++)
        {
            this.balls[i].transform.position = new Vector3(this.pos[i * 3], this.pos[i * 3 + 1], this.pos[i * 3 + 2]);
        }
    }

    public void simulate(float dt, Vector3 gravity, float[] worldBounds)
    {
        var minDist = 2.0f * this.radius;

        // 重力影响
        for (var i = 0; i < this.numBalls; i++)
        {
            vecAdd(this.vel, i, new float[] { gravity.x, gravity.y, gravity.z}, 0, dt);
            vecCopy(this.prevPos, i, this.pos, i);
            vecAdd(this.pos, i, this.vel, i, dt);
        }
        this.hash.create(this.pos);

        // 处理碰撞
        for (var i = 0; i < this.numBalls; i++)
        {
            this.balls[i].GetComponent<MeshRenderer>().material.color = Color.green;
            // 墙壁碰撞
            for (var dim = 0; dim < 3; dim++)
            {

                var nr = 3 * i + dim;
                if (this.pos[nr] < worldBounds[dim] + this.radius)
                {
                    this.pos[nr] = worldBounds[dim] + this.radius;
                    this.vel[nr] = -this.vel[nr];
                    this.balls[i].GetComponent<MeshRenderer>().material.color = Color.yellow;
                }
                else if (this.pos[nr] > worldBounds[dim + 3] - this.radius)
                {
                    this.pos[nr] = worldBounds[dim + 3] - this.radius;
                    this.vel[nr] = -this.vel[nr];
                    this.balls[i].GetComponent<MeshRenderer>().material.color = Color.yellow;
                }
            }

            //  球体相互碰撞
            this.hash.query(this.pos, i, 2.0f * this.radius);

            for (var nr = 0; nr < this.hash.querySize; nr++)
            {
                var j = this.hash.queryIds[nr];

                vecSetDiff(this.normal, 0, this.pos, i, this.pos, j);
                var d2 = vecLengthSquared(this.normal, 0);

                if (d2 > 0.0 && d2 < minDist * minDist)
                {
                    float d = Mathf.Sqrt(d2);
                    vecScale(this.normal, 0, 1.0f / d);

                    var corr = (minDist - d) * 0.5f;
                    vecAdd(this.pos, i, this.normal, 0, corr);
                    vecAdd(this.pos, j, this.normal, 0, -corr);


                    var vi = vecDot(this.vel, i, this.normal, 0);
                    var vj = vecDot(this.vel, j, this.normal, 0);

                    vecAdd(this.vel, i, this.normal, 0, vj - vi);
                    vecAdd(this.vel, j, this.normal, 0, vi - vj);

                    this.balls[i].GetComponent<MeshRenderer>().material.color = Color.yellow;
                }
            }
        }
        this.updateMesh();
    }
}

public class lesson_ele : MonoBehaviour
{
    ThreeDBalls balls;
    float[] s = { 0, 0, 0, 1, 1, 1 };

    public void initPhysics()
    {
        var radius = 0.025f;

        var spacing = 3.0f * radius;
        var velRand = 0.2f;

        var numX = Mathf.FloorToInt(Mathf.Floor((s[3] - s[0] - 2.0f * spacing) / spacing));
        var numY = Mathf.FloorToInt(Mathf.Floor((s[4] - s[1] - 2.0f * spacing) / spacing));
        var numZ = Mathf.FloorToInt(Mathf.Floor((s[5] - s[2] - 2.0f * spacing) / spacing));

        var pos = new float[3 * numX * numY * numZ];
        var vel = new float[3 * numX * numY * numZ];
        //vel.fill(0.0);

        for (var xi = 0; xi < numX; xi++)
        {
            for (var yi = 0; yi < numY; yi++)
            {
                for (var zi = 0; zi < numZ; zi++)
                {
                    var x = 3 * ((xi * numY + yi) * numZ + zi);
                    var y = x + 1;
                    var z = x + 2;
                    pos[x] = s[0] + spacing + xi * spacing;
                    pos[y] = s[1] + spacing + yi * spacing;
                    pos[z] = s[2] + spacing + zi * spacing;

                    vel[x] = -velRand + 2.0f * velRand * Random.Range(0f, 1f);
                    vel[y] = -velRand + 2.0f * velRand * Random.Range(0f, 1f);
                    vel[z] = -velRand + 2.0f * velRand * Random.Range(0f, 1f);
                }
            }
        }

        balls = new ThreeDBalls(radius, pos, vel);

    }
    void Start()
    {
        initPhysics();
    }

    private void FixedUpdate()
    {
        balls.simulate(Time.fixedDeltaTime, Vector3.down * .98f, s);
    }
}
