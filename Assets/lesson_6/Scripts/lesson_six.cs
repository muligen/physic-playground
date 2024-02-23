using PhysicShapes;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Pendulum
{
    public List<float> masses;
    public List<float> lengths;
    public List<Vector2> pos;
    public List<Vector2> prevPos;
    public List<Vector2> vel;
    public Pendulum(List<float> m, List<float> l, List<Vector2> p)
    {
        this.masses = m;
        this.lengths = l;
        this.pos = p;
        this.prevPos = new List<Vector2>();
        this.vel = new List<Vector2>();
        for (int i = 0; i < masses.Count; i++)
        {
            this.prevPos.Add(Vector2.zero);
            this.vel.Add(Vector2.zero);
        }
    }
    public void simulate(float dt, Vector2 gravity){
        var p = this;
        for (var i = 1; i < p.masses.Count; i++)
        {
            p.vel[i] += dt * gravity;
            p.prevPos[i] = p.pos[i];
            p.pos[i] += p.vel[i] * dt;
        }
        for (var i = 1; i < p.masses.Count; i++)
        {
            Vector2 dir = p.pos[i] - p.pos[i - 1];
            float d = dir.magnitude;
            dir = dir / d;
            float w0 = p.masses[i - 1] > 0.0f ? 1.0f / p.masses[i - 1] : 0.0f;
            float w1 = p.masses[i] > 0.0f ? 1.0f / p.masses[i] : 0.0f;
            float corr = (p.lengths[i - 1] - d) / (w0 + w1);
            p.pos[i - 1] -= w0 * dir * corr;
            p.pos[i] += w1 * dir * corr;
        }
        for (var i = 1; i < p.masses.Count; i++)
        {
            p.vel[i] = (p.pos[i] - p.prevPos[i]) / dt;
        }

    }
}


public class lesson_six : MonoBehaviour
{
    public Image img_ball1;
    public Image img_ball2;
    public Image img_ball3;
    public Image img_ball4;

    public float mass1 = 0f;
    public float mass2 = 100f;
    public float mass3 = 50f;
    public float mass4 = 100f;

    private LineRenderer lineRenderer;

    Pendulum p = null;
    
    List<Bead> list_bead;
    void Start()
    {
        p = new Pendulum(
            new List<float>() { mass1, mass2, mass3, mass4 },
            new List<float>() { 200, 200, 200 },
            new List<Vector2>() { img_ball1.transform.position, img_ball2.transform.position, img_ball3.transform.position, img_ball4.transform.position }
        );

        //Ìí¼ÓLineRenderer×é¼þ
        lineRenderer = gameObject.AddComponent<LineRenderer>();
        lineRenderer.positionCount = 4;
        lineRenderer.startWidth = 5f;
    } 
    private void FixedUpdate()
    {

        float dt = Time.fixedDeltaTime;

        var sdt = dt / 100;

        for (var step = 0; step < 100; step++)
        {
            p.simulate(sdt, Vector2.down * 9.8f * 200);
        }
        img_ball1.transform.position = p.pos[0];
        img_ball2.transform.position = p.pos[1];
        img_ball3.transform.position = p.pos[2];
        img_ball4.transform.position = p.pos[3];

        img_ball3.rectTransform.sizeDelta = img_ball2.rectTransform.sizeDelta * MathF.Sqrt(p.masses[2] / p.masses[1]);
        img_ball4.rectTransform.sizeDelta = img_ball2.rectTransform.sizeDelta * MathF.Sqrt(p.masses[3] / p.masses[1]);

        lineRenderer.SetPosition(0, new Vector3(p.pos[0].x, p.pos[0].y, 0));
        lineRenderer.SetPosition(1, new Vector3(p.pos[1].x, p.pos[1].y, 0));
        lineRenderer.SetPosition(2, new Vector3(p.pos[2].x, p.pos[2].y, 0));
        lineRenderer.SetPosition(3, new Vector3(p.pos[3].x, p.pos[3].y, 0));
    }
}
