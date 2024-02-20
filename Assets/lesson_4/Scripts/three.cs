using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using UnityEngine.UI;
using PhysicShapes;


public class Flipper {
    public float radius;
    public Vector2 pos;
    public float length;
    public float restAngle;
    public float maxRotation;
    public int sign;
    public float angularVelocity;
    public float rotation;
    public float currentAngularVelocity;
    public Flipper(float radius, Vector2 pos, float length, float restAngle, float maxRotation, float angularVelocity)
    {
        this.radius = radius;
        this.pos = pos;
        this.length = length;
        this.restAngle = restAngle;
        this.maxRotation = Mathf.Abs(maxRotation);
        this.sign = Math.Sign(maxRotation);
        this.angularVelocity = angularVelocity;
        // changing
        this.rotation = 0;
        this.currentAngularVelocity = 0;
        //this.touchIdentifier = -1;
    }

    public void simulate(float dt)
    {
        var prevRotation = this.rotation;
        if (Input.GetMouseButton(0) && select(Input.mousePosition))
            this.rotation = Math.Min(this.rotation + dt * this.angularVelocity, this.maxRotation);
        else
            this.rotation = Math.Max(this.rotation - dt * this.angularVelocity, 0.0f);
        this.currentAngularVelocity = this.sign * (this.rotation - prevRotation) / dt;
    }
    public bool select(Vector3 pos)
    {
        var d = new Vector2(pos.x, pos.y);
        d = this.pos - d;
        return d.magnitude < this.length + 50f;
    }

    public Vector2 getTip()
    {
        var angle = this.restAngle + this.sign * this.rotation;
        var dir = new Vector2((float)Math.Cos(angle), (float)Math.Sin(angle));
        var tip = this.pos + this.length * dir;
        return tip;
    }

    public float getDegree()
    {
        var angle = this.restAngle + this.sign * this.rotation;
        return angle * Mathf.Rad2Deg;
    }
}

public class Obstacle
{
    public Vector2 pos;
    public float radius;
    public float pushVel;
    public Obstacle(Vector2 pos, float radius, float pushVel)
    {
        this.pos = pos;
        this.radius = radius;
        this.pushVel = pushVel;
    }
}



public class three : MonoBehaviour
{
    public Image ImgBall1;
    public Image ImgBall2;
    TwoDBall ball1;
    TwoDBall ball2;
    public Image ImgF1;
    public Image ImgF2;
    Flipper f1;
    Flipper f2;
    public Image ImgObstacle1;
    public Image ImgObstacle2;
    public Image ImgObstacle3;
    Obstacle obstacle1;
    Obstacle obstacle2;
    Obstacle obstacle3;
    public Text score;

    List<TwoDBall> BallList;
    List<Vector2> borders;
    // Start is called before the first frame update
    void Start()
    {
        ball1 = new TwoDBall(Vector2.zero, 25f, ImgBall1.transform.position, 10, 0.5f);
        ball2 = new TwoDBall(Vector2.zero, 25f, ImgBall2.transform.position, 10, 0.5f);
        f1 = new Flipper(30f, ImgF1.transform.position, 180f, ImgF1.transform.rotation.eulerAngles.z * Mathf.Deg2Rad, 1, 10.0f);
        f2 = new Flipper(30f, ImgF2.transform.position, 180f, ImgF2.transform.rotation.eulerAngles.z * Mathf.Deg2Rad, -1, 10.0f);
        obstacle1 = new Obstacle(ImgObstacle1.transform.position, 50f, 1000f);
        obstacle2 = new Obstacle(ImgObstacle2.transform.position, 50f, 1000f);
        obstacle3 = new Obstacle(ImgObstacle3.transform.position, 50f, 1000f);
        BallList = new List<TwoDBall>() { ball1, ball2};

        //初始化墙壁
        borders = new List<Vector2>()
        {
            new Vector2(500, 100),
            new Vector2(800, 100),
            new Vector2(800, 300),
            new Vector2(1000, 350),
            new Vector2(1000, 1050),
            new Vector2(300, 1050),
            new Vector2(300, 350),
            new Vector2(500, 300),
        };
    }

    private void Update()
    {
        //画墙
        for (int i = 0; i < borders.Count; i++)
        {
            Debug.DrawLine(borders[i], borders[(i + 1) % borders.Count]);
        }
    }

    private void FixedUpdate()
    {
        f1.simulate(Time.fixedDeltaTime);
        f2.simulate(Time.fixedDeltaTime);
        for (int i = 0; i < BallList.Count; i++)
        {
            BallList[i].Simulate(9.8f * 200, Time.fixedDeltaTime);
            for (int j = i+1; j < BallList.Count; j++)
            {
                BallBallCollsion(BallList[i], BallList[j]);
            }
            BallFlipperCollsion(BallList[i], f1);
            BallFlipperCollsion(BallList[i], f2);
            BallObstacleCollsion(BallList[i], obstacle1);
            BallObstacleCollsion(BallList[i], obstacle2);
            BallObstacleCollsion(BallList[i], obstacle3);
            BallWallCollsion(BallList[i], borders);
        }

        ImgF1.transform.rotation = Quaternion.Euler(0, 0, f1.getDegree());
        ImgF2.transform.rotation = Quaternion.Euler(0, 0, f2.getDegree());
        ImgBall1.transform.position = ball1.pos;
        ImgBall2.transform.position = ball2.pos;
    }


    public Vector2 closestPointOnSegment(Vector2 p, Vector2 a, Vector2 b) //寻找点到胶囊体最近点
    {
        var ab = b - a;
        var k = Vector2.Dot(ab, ab);
        if (k == 0.0)
            return a;
        var t = (Vector2.Dot(p, ab) - Vector2.Dot(a, ab)) / k;
        t = Math.Max(0, Math.Min(1, t));
        return a + t * ab;
    }
    public void BallBallCollsion(TwoDBall ball1, TwoDBall ball2)
    {
        var restitution = Math.Min(ball1.restitution, ball2.restitution);
        var dir = new Vector2();
        dir = ball2.pos - ball1.pos;
        var d = dir.magnitude;
        if (d == 0.0 || d > ball1.radius + ball2.radius) //去掉除0问题
            return;

        dir /= d;

        float corr = (ball1.radius + ball2.radius - d) / 2;
        ball1.pos -= corr * dir;
        ball2.pos += corr * dir;

        float v1 = Vector2.Dot(ball1.vel, dir);
        float v2 = Vector2.Dot(ball2.vel, dir);

        var m1 = ball1.mass;
        var m2 = ball2.mass;

        var newV1 = (m1 * v1 + m2 * v2 - m2 * (v1 - v2) * restitution) / (m1 + m2);
        var newV2 = (m1 * v1 + m2 * v2 - m1 * (v2 - v1) * restitution) / (m1 + m2); 

        ball1.vel += (newV1 - v1) * dir; //碰撞方向速度归0所以去掉v1，在赋予新速度newV1。
        ball2.vel += (newV2 - v2) * dir; 
    }
    public void BallObstacleCollsion(TwoDBall ball, Obstacle obstacle)
    {
        var dir = new Vector2();
        dir = ball.pos - obstacle.pos;
        var d = dir.magnitude;
        if (d == 0.0 || d > ball.radius + obstacle.radius)
            return;

        dir /= d;

        var corr = ball.radius + obstacle.radius - d; //这里没有添加回弹距离？
        ball.pos += corr * dir;

        var v = Vector2.Dot(ball.vel, dir);
        ball.vel += (obstacle.pushVel - v) * dir;

        score.text = (int.Parse(score.text) + 1).ToString();
    }

    public void BallWallCollsion(TwoDBall ball, List<Vector2> border)
    {
        if (border.Count < 3)
            return;

        // find closest segment;
        var d = new Vector2();
        var closest = new Vector2();
        var ab = new Vector2();
        var normal = new Vector2();
        var minDist = 0.0;
        for (var i = 0; i < border.Count; i++)
        {
            var a = border[i];
            var b = border[(i + 1) % border.Count];
            var c = closestPointOnSegment(ball.pos, a, b);
            d = ball.pos - c;
            if (i == 0 || d.magnitude < minDist)
            {
                minDist = d.magnitude;
                closest = c;
                ab = b - a;
                normal = perp(ab);
            }
        }

        // push out
        d = ball.pos - closest;
        var dist = d.magnitude;
        if (dist == 0.0)
        {
            d = normal;
            dist = normal.magnitude;
        }
        d /= dist;
        if (Vector2.Dot(d, normal) >= 0.0)
        {
            if (dist > ball.radius)
                return;
            ball.pos += (ball.radius - dist) * d;
        }
        else
            ball.pos -= (dist + ball.radius) * d;
        // update velocity
        var v = Vector2.Dot(d, ball.vel);
        var vnew = Math.Abs(v) * ball.restitution;

        ball.vel += (vnew - v) * d;
    }

    public Vector2 perp(Vector2 v)
    {
        return new Vector2(-v.y, v.x);
    }
    public void BallFlipperCollsion(TwoDBall ball, Flipper flipper)
    {
        var closest = closestPointOnSegment(ball.pos, flipper.pos, flipper.getTip());
        var dir = ball.pos - closest;
        var d = dir.magnitude;
        if (d == 0.0 || d > ball.radius + flipper.radius)
            return;

        dir = dir.normalized;

        var corr = (ball.radius + flipper.radius - d);
        ball.pos += corr * dir;

        //update velocitiy

        var radius = closest;
        radius += flipper.radius * dir;
        radius -= flipper.pos;
        var surfaceVel = perp(radius);
        surfaceVel *= flipper.currentAngularVelocity;

        var v = Vector2.Dot(ball.vel, dir);
        var vnew = Vector2.Dot(surfaceVel, dir);

        ball.vel += (vnew - v) * dir;
    }


}
