import rclpy
from rclpy.node import Node
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint

class TrajectoryPublisher(Node):
    def __init__(self):
        super().__init__('trajectory_publisher')
        self.publisher = self.create_publisher(JointTrajectory, '/joint_trajectory_controller/joint_trajectory', 10)
        self.timer = self.create_timer(1.0, self.publish_trajectory)

    def publish_trajectory(self):
        msg = JointTrajectory()
        msg.joint_names = ['Revolute9', 'revolute10', 'Revolute11']  # List of joint names
        point = JointTrajectoryPoint()
        point.positions = [1.57, 0.57, 1.0]  # List of joint positions
        point.time_from_start = rclpy.duration.Duration(seconds=10)  # Time duration
        msg.points.append(point)
        self.publisher.publish(msg)

def main(args=None):
    rclpy.init(args=args)
    trajectory_publisher = TrajectoryPublisher()
    rclpy.spin(trajectory_publisher)
    trajectory_publisher.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()

