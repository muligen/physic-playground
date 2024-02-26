using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class SoftBody
{
    public int numParticles;
    public int numTets;
    public float[] pos;
    public float[] prevPos;
    public float[] vel;

    public int[] tetIds;
    public int[] edgeIds;

    public float[] restVol;
    public float[] edgeLengths;
    public float[] invMass;

    public float edgeCompliance;
    public float volCompliance;

    public float[] temp;
    public float[] grads;

    public int grabId;
    public float grabInvMass;

    public GameObject surfaceMesh;
    public int[,] volIdOrder = { { 1, 3, 2 }, { 0, 2, 3 }, { 0, 3, 1 }, { 0, 1, 2 } };


    public SoftBody(BunnyDataStruct tetMesh, int scene, float edgeCompliance = 100.0f, float volCompliance = 0.0f)
    {
        this.numParticles = tetMesh.verts.Count / 3;
        this.numTets = tetMesh.tetIds.Count / 4;
        this.pos = tetMesh.verts.ToArray(); // 拷贝构造
        this.prevPos = tetMesh.verts.ToArray(); // 拷贝构造
        this.vel = new float[3 * this.numParticles];

        this.tetIds = tetMesh.tetIds.ToArray();
        this.edgeIds = tetMesh.tetEdgeIds.ToArray();
        this.restVol = new float[this.numTets]; // 定义固定长度的List，count = 0？
        this.edgeLengths = new float[this.edgeIds.Length / 2];
        this.invMass = new float[this.numParticles];

        this.edgeCompliance = edgeCompliance;
        this.volCompliance = volCompliance;

        this.temp = new float[4 * 3];
        this.grads = new float[4 * 3];

        this.grabId = -1;
        this.grabInvMass = 0.0f;

        this.initPhysics(); 

        // surface tri mesh  三角面渲染。
        surfaceMesh = new GameObject("Cube");
        surfaceMesh.transform.position = Vector3.zero;

        //顶点数组
        Vector3[] _vertices = new Vector3[this.numParticles];
        for (int i = 0; i < this.numParticles; i++)
        {
            _vertices[i] = new Vector3(this.pos[i * 3], this.pos[i * 3 + 1], this.pos[i * 3 + 2]);
        }
        int[] _triangles = tetMesh.tetSurfaceTriIds.ToArray();
        Mesh mesh = new Mesh()
         {
            vertices = _vertices,
            triangles = _triangles,
         };
        mesh.RecalculateNormals();
        surfaceMesh.AddComponent<MeshFilter>().mesh = mesh;
        var a = Resources.Load<Material>("bunny/mat_pink");
        surfaceMesh.AddComponent<MeshRenderer>().material = a;
    }

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
    public void translate(float x, float y, float z)
    {
        for (var i = 0; i < this.numParticles; i++)
        {
            var t = new float[] { x, y, z };
            vecAdd(this.pos, i, t, 0);
            vecAdd(this.prevPos, i, t, 0);
        }
    }

    public void updateMeshes()
    {
        var mesh = this.surfaceMesh.GetComponent<MeshFilter>().mesh;
        Vector3[] _vertices = new Vector3[this.numParticles];
        for (int i = 0; i < this.numParticles; i++)
        {
            _vertices[i] = new Vector3(this.pos[i * 3], this.pos[i * 3 + 1], this.pos[i * 3 + 2]);
        }
        mesh.vertices = _vertices;
        mesh.RecalculateNormals();
    }

    public float getTetVolume(int nr)
    {
        var id0 = this.tetIds[4 * nr];
        var id1 = this.tetIds[4 * nr + 1];
        var id2 = this.tetIds[4 * nr + 2];
        var id3 = this.tetIds[4 * nr + 3];
        vecSetDiff(this.temp, 0, this.pos, id1, this.pos, id0);
        vecSetDiff(this.temp, 1, this.pos, id2, this.pos, id0);
        vecSetDiff(this.temp, 2, this.pos, id3, this.pos, id0);
        vecSetCross(this.temp, 3, this.temp, 0, this.temp, 1);
        return vecDot(this.temp, 3, this.temp, 2) / 6.0f;
    }

    public void initPhysics()
    {
        for (int i = 0; i < this.invMass.Length; i++)
        {
            this.invMass[i] = 0;
        }
        for (int i = 0; i < this.restVol.Length; i++) { 
            this.restVol[i] = 0;
        }
        // 每一个体积，平均分摊质量mass给4个顶点。
        for (var i = 0; i < this.numTets; i++)
        {
            var vol = this.getTetVolume(i);
            this.restVol[i] = vol;
            var pInvMass = vol > 0.0f ? 1.0f / (vol / 4.0f) : 0.0f;
            this.invMass[this.tetIds[4 * i]] += pInvMass;
            this.invMass[this.tetIds[4 * i + 1]] += pInvMass;
            this.invMass[this.tetIds[4 * i + 2]] += pInvMass;
            this.invMass[this.tetIds[4 * i + 3]] += pInvMass;
        }
        // 计算边长
        for (var i = 0; i < this.edgeLengths.Length; i++)
        {
            var id0 = this.edgeIds[2 * i];
            var id1 = this.edgeIds[2 * i + 1];
            this.edgeLengths[i] = Mathf.Sqrt(vecDistSquared(this.pos, id0, this.pos, id1));
        }
    }

    public void preSolve(float dt, Vector3 gravity)
    {
        for (var i = 0; i < this.numParticles; i++)
        {
            if (this.invMass[i] == 0.0)
                continue;
            vecAdd(this.vel, i, new float[] { gravity.x, gravity.y, gravity.z}, 0, dt);
            vecCopy(this.prevPos, i, this.pos, i);
            vecAdd(this.pos, i, this.vel, i, dt);
            var y = this.pos[3 * i + 1]; // 限制y分量不低于0，地面
            if (y < 0.0)
            {
                vecCopy(this.pos, i, this.prevPos, i);
                this.pos[3 * i + 1] = 0.0f;
            }
        }
    }

    public void solve(float dt)
    {
        this.solveEdges(this.edgeCompliance, dt);
        this.solveVolumes(this.volCompliance, dt);
    }

    public void postSolve(float dt)
    {
        for (var i = 0; i < this.numParticles; i++)
        {
            if (this.invMass[i] == 0.0)
                continue;
            vecSetDiff(this.vel, i, this.pos, i, this.prevPos, i, 1.0f / dt);
        }
        this.updateMeshes();
    }

    public void solveEdges(float compliance, float dt)
    {
        var alpha = compliance / dt / dt;  // xpbd的区别，加了这个项

        for (var i = 0; i < this.edgeLengths.Length; i++)
        {
            var id0 = this.edgeIds[2 * i];
            var id1 = this.edgeIds[2 * i + 1];
            var w0 = this.invMass[id0];
            var w1 = this.invMass[id1];
            var w = w0 + w1;
            if (w == 0.0)
                continue;

            vecSetDiff(this.grads, 0, this.pos, id0, this.pos, id1);
            var len = Mathf.Sqrt(vecLengthSquared(this.grads, 0));
            if (len == 0.0)
                continue;
            vecScale(this.grads, 0, 1.0f / len);
            var restLen = this.edgeLengths[i];
            var C = len - restLen;
            var s = -C / (w + alpha);
            vecAdd(this.pos, id0, this.grads, 0, s * w0);
            vecAdd(this.pos, id1, this.grads, 0, -s * w1);
        }
    }

    public void solveVolumes(float compliance, float dt)
    {
        var alpha = compliance / dt / dt;

        for (var i = 0; i < this.numTets; i++)
        {
            var w = 0.0f;

            for (var j = 0; j < 4; j++)
            {
                var id0 = this.tetIds[4 * i + this.volIdOrder[j, 0]];
                var id1 = this.tetIds[4 * i + this.volIdOrder[j, 1]];
                var id2 = this.tetIds[4 * i + this.volIdOrder[j, 2]];

                vecSetDiff(this.temp, 0, this.pos, id1, this.pos, id0);
                vecSetDiff(this.temp, 1, this.pos, id2, this.pos, id0);
                vecSetCross(this.grads, j, this.temp, 0, this.temp, 1);
                vecScale(this.grads, j, 1.0f / 6.0f);

                w += this.invMass[this.tetIds[4 * i + j]] * vecLengthSquared(this.grads, j);
            }
            if (w == 0.0)
                continue;

            var vol = this.getTetVolume(i);
            var restVol = this.restVol[i];
            var C = vol - restVol;
            var s = -C / (w + alpha);

            for (var j = 0; j < 4; j++)
            {
                var id = this.tetIds[4 * i + j];
                vecAdd(this.pos, id, this.grads, j, s * this.invMass[id]);
            }
        }
    }

    public void squash()
    {
        for (var i = 0; i < this.numParticles; i++)
        {
            this.pos[3 * i + 1] = 0.5f;
        }
        this.updateMeshes();
    }

    public void startGrab(Vector3 pos)
    {
        float[] p = { pos.x, pos.y, pos.z };
        var minD2 = float.MaxValue;
        this.grabId = -1;
        for (int i = 0; i < this.numParticles; i++)
        {
            var d2 = vecDistSquared(p, 0, this.pos, i);
            if (d2 < minD2)
            {
                minD2 = d2;
                this.grabId = i;
            }
        }

        if (this.grabId >= 0)
        {
            this.grabInvMass = this.invMass[this.grabId];
            this.invMass[this.grabId] = 0.0f;
            vecCopy(this.pos, this.grabId, p, 0);
        }
    }

    public void moveGrabbed(Vector3 pos, Vector3 vel)
    {
        if (this.grabId >= 0)
        {
            float[] p = { pos.x, pos.y, pos.z };
            vecCopy(this.pos, this.grabId, p, 0);
        }
    }

    public void endGrab(Vector3 pos, Vector3 vel)
    {
        if (this.grabId >= 0)
        {
            this.invMass[this.grabId] = this.grabInvMass;
            float[] v = { vel.x, vel.y, vel.z };
            vecCopy(this.vel, this.grabId, v, 0);
        }
        this.grabId = -1;
    }
}

public class lesson_ten : MonoBehaviour
{
    List<SoftBody> softBodies;
    public Slider slider;
    public Button btn;

    public void OnClick_btn()
    {
        var d = new BunnyData();
        var bunny = new SoftBody(d.data, 0);
        softBodies.Add(bunny);
    }
    void Awake()
    {
        var d = new BunnyData();
        var bunny = new SoftBody(d.data, 0);
        softBodies = new List<SoftBody>
        {
            bunny
        };
        btn.onClick.AddListener(OnClick_btn);
    }
    void simulate(float dt)
    {
        int numSubsteps = 10; //substep太少导致模拟崩坏，错误。
        var sdt = dt / numSubsteps;

        for (var step = 0; step < numSubsteps; step++)
        {

            for (var i = 0; i < softBodies.Count; i++)
                softBodies[i].preSolve(sdt, Vector3.down * 10f);

            for (var i = 0; i < softBodies.Count; i++)
                softBodies[i].solve(sdt);

            for (var i = 0; i < softBodies.Count; i++)
                softBodies[i].postSolve(sdt);

        }
        if (Input.GetMouseButtonDown(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            RaycastHit hit;
            if (Physics.Raycast(ray, out hit, Mathf.Infinity))
            {
                softBodies[0].startGrab(hit.point);
            }
        }
        if (Input.GetMouseButton(0) )
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            RaycastHit hit;
            if (Physics.Raycast(ray, out hit, Mathf.Infinity))
            {
                //Debug.Log("坐标: " + hit.point);
                softBodies[0].moveGrabbed(hit.point, Vector3.zero);
            }
        }
        else
        {
            softBodies[0].endGrab(Vector3.zero, Vector3.zero);
        }

        softBodies[0].edgeCompliance = slider.value * 50f;
    }
    private void FixedUpdate()
    {
        simulate(Time.fixedDeltaTime);
    }
}
